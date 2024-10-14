class PRINTS {
    static parseOutput(
            String printsOutput,
            String hierarchyDb
    ) {
        // Build up a map of the fingerprint hierarchies
        Map<String, HierarchyEntry> hierarchyMap = new LinkedHashMap<>()
        File hierarchyFile = new File(hierarchyDb)
        hierarchyFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (!line.startsWith("#")) {
                    def row = line.trim().split("\\|")
                    if (row.length >= 3) {
                        // e.g. G6PDHDRGNASE|PR00079|1e-04|0|
                        final String modelId = row[0].trim()
                        final String modelAccession = row[1].trim()
                        final Double evalueCutoff = Double.parseDouble(row[2].trim())
                        final int minMotifCount = Integer.parseInt(row[3].trim())
                        HierarchyEntry hierarchyEntry = new HierarchyEntry(modelId, modelAccession, evalueCutoff, minMotifCount)
                        if (row.length > 4) {
                            // e.g. DHBDHDRGNASE|PR01397|1e-04|0|SDRFAMILY or TYROSINASE|PR00092|1e-04|0|*
                            String siblingsOrDomain = row[4].trim()
                            if (siblingsOrDomain == "*") {
                                boolean isDomain = true
                                hierarchyEntry.changeDomainStatus(isDomain)
                            } else if (siblingsOrDomain.size() > 0) {
                                String[] siblingsIds = siblingsOrDomain.split("\\,")
                                hierarchyEntry.addSiblings(siblingsIds)
                            }
                        }
                        hierarchyMap.put(modelId, hierarchyEntry)
                    }
                }
            }
        }

        // parse the Prints Output into Match instances
        File printsFile = new File(printsOutput)
        if (!printsFile.exists()){
            System.out.println("Could not find Hierarchy DB for PRINTS")
            System.exit 1
        }

        Map<String, Map<String, Match>> hits = new LinkedHashMap<>()  // <protein ID <Model Acc, Match>>
        Map<String, String> name2accession = new LinkedHashMap<>()  // <motifName, motifAccession>
        String queryAccession = null  // protein seq ID
        Map<String, Match> thisProteinsMatches = new LinkedHashMap<>()  // all Matches for the current working protein

        printsFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Sn;")) {
                    // Post-process/filter the hits from the last protein
                    if (!thisProteinsMatches.isEmpty()) {
                        hits = filterProteinMatches(queryAccession, thisProteinsMatches, hierarchyMap, name2accession, hits)
                    }
                    // Start the new protein: Get the query sequence id
                    queryAccession = line.replaceFirst("Sn; ", "").trim()
                    name2accession = new LinkedHashMap<>()
                    thisProteinsMatches = new LinkedHashMap<>()
                }
                // Do not retrieve the motif Description from the 1TBH line as this is retrieved in the XREFS subworkflow
                // And retrieve the motifAcc from the heirarchyDB instead of the 1TBH line
                else if (line.startsWith("2TBH")) {  // Used to create Match instances
                    // Line: 2TBH  modelId  NumMotifs  SumId  AveId  ProfScore  Ppvalue  Evalue  GraphScan
                    def matcher = line =~ ~/^2TBH\s+(.+?)\s+(\d)\s+of\s+\d\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\s+.*?\s+(.*?)\s+([\.iI]+)\s+$/
                    if (matcher.find()) {
                        String modelId = matcher.group(1).trim()
                        int numOfMotifs = Integer.parseInt(matcher.group(2).trim())
                        Double evalue = Double.parseDouble(matcher.group(3).trim())
                        String graphScan = matcher.group(4).trim()

                        if (hierarchyMap.containsKey(modelId)) {
                            String modelAccession = hierarchyMap[modelId]["modelAccession"]
                            int minMotifCount = hierarchyMap[modelId]["minMotifCount"] as int
                            Double cutoff = hierarchyMap[modelId]["evalueCutoff"] as Double

                            if (evalue <= cutoff && numOfMotifs > minMotifCount) {
                                thisProteinsMatches.put(modelAccession, new Match(modelAccession, evalue, graphScan))
                                name2accession.put(modelId, modelAccession)
                            }
                        }
                    }
                }
                else if (line.startsWith("3TBH")) {  // Used to create Location instances
                    // Line: 3TBH modelId NoOfMotifs IdScore PfScore Pvalue Sequence Len Low Pos High
                    def matcher = line =~ ~/^3TBH\s+(.+?)\s+(\d+)\s+of\s+\d+\s+([\d\.]+)\s+\d+\s+(.+?)\s+([\w#]+)\s+(\d+)\s+\d+\s+(-?\d+)\s*\d\s*$/
                    if (matcher.find()) {
                        // Check motif passed the e-value cutoff when parsing line 2TBH - else SKIP!
                        String modelId = matcher.group(1).trim()
                        String modelAccession = name2accession[modelId]
                        if (thisProteinsMatches.containsKey(modelAccession)) {
                            final int motifNumber = Integer.parseInt(matcher.group(2).trim())
                            final Float idScore = Float.parseFloat(matcher.group(3).trim())
                            final Double pvalue = Double.parseDouble(matcher.group(4).trim())
                            final String sequence = matcher.group(5).trim()
                            int motifLength = Integer.parseInt(matcher.group(6).trim())
                            int position = Integer.parseInt(matcher.group(7).trim())

                            int end = position + motifLength - 1
                            if (position < 1) {
                                position = 1
                            }

                            // If the motif seq ends in #, it over hands the end of the query sequence
                            // If that's the case, find out by how far and adjust the end appropriately
                            int indexCheck = 0
                            if (sequence.endsWith("#")) {
                                motifLength = sequence.length()
                                indexCheck = motifLength - 1
                                while (sequence[indexCheck] == "#") {
                                    indexCheck -= 1
                                }
                            }
                            end = end - (motifLength - indexCheck) + 1

                            Location location = new Location(position, end, pvalue, idScore, motifNumber)
                            thisProteinsMatches[modelAccession].addLocation(location)

                        }
                    }
                }
            }
        }
        // process the last protein
        if (!thisProteinsMatches.isEmpty()) {
            hits = filterProteinMatches(queryAccession, thisProteinsMatches, hierarchyMap, name2accession, hits)
        }

        return hits
    }

    static filterProteinMatches(
            String queryAccession,
            Map<String, Match> thisProteinMatches,
            Map<String, HierarchyEntry> hierarchyMap,
            Map<String, String> name2accession,
            Map<String, Map<String, Match>> hits
    ) {
        List<Match> sortedMatches = sortMatches(thisProteinMatches.values() as List<Match>)
        Set<String> HierarchtEntryIdLimitation = hierarchyMap.keySet() // initialise with all model IDs
        for (Match match in sortedMatches) {
            String modelId = name2accession.find{ it.value == match.modelAccession } ?.key
            HierarchyEntry hierarchyEntry = hierarchyMap[modelId]
            // check if match passes the filtering criteria
            if (selectMatch(match, hierarchyEntry, HierarchtEntryIdLimitation)) {
                if (!hits.containsKey(queryAccession)) {hits.put(queryAccession, new LinkedHashMap<>())}
                hits[queryAccession].put(match.modelAccession, match)
                // passed the filter and may have its own hierarchy so update the hierarchy limitation
                HierarchtEntryIdLimitation = hierarchyEntry.siblingsIds
            }
        }

        return hits
    }

    static sortMatches(List<Match> matches) { // This comparator is CRITICAL to the working of PRINTS post-processing
        return matches.sort { matchA, matchB ->
            int evalueComparison = matchA.evalue <=> matchB.evalue
            if (evalueComparison != 0) return evalueComparison

            int modelAccessionComparison = matchA.modelAccession <=> matchB.modelAccession
            if (modelAccessionComparison != 0) return modelAccessionComparison
            // NOTE: There is only one location per match
            int motifNumberComparison = matchA.locations[0].motifNumber <=> matchB.locations[0].motifNumber
            if (motifNumberComparison != 0) return motifNumberComparison

            int startLocationComparison = matchA.locations[0].start <=> matchB.locations[0].start
            if (startLocationComparison != 0) return startLocationComparison

            return matchA.locations[0].end <=> matchB.locations[0].end
        }
    }

    static selectMatch(Match match, HierarchyEntry hierarchy, Set<String> HierarchyEntryIdLimitation) {
        // if a domain: PASS - there is no hierarchy to deal with
        if (hierarchy.isDomain) {return true}
        // if the previous model limited the filtering to its siblings. If the first model, this is all models in HierarchyDB
        if (HierarchyEntryIdLimitation.contains(match.modelAccession)) {return true}
    }

}

class HierarchyEntry implements Serializable { // represents a row in the HierarchyDB
    String modelId
    String modelAccession
    Double evalueCutoff
    int minMotifCount
    boolean isDomain = false
    String[] siblingsIds = null
    
    HierarchyEntry(String modelId, String modelAccession, Double evalueCutoff, int minMotifCount) {
        this.modelId = modelId
        this.modelAccession = modelAccession
        this.evalueCutoff = evalueCutoff
        this.minMotifCount = minMotifCount
    }

    void changeDomainStatus(boolean isDomain) {
        this.isDomain = isDomain
    }

    void addSiblings(String[] siblingIds) {
        this.siblingsIds = siblingsIds
    }
}
