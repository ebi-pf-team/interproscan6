class PRINTS {
    static parseOutput(
            String printsOutput,
            String hierarchyDb
    ) {
        // Build up a map of the Model ID to fingerprint hierarchies
        Map<String, HierarchyEntry> hierarchyMap = HierarchyEntry.parseHierarchyDbFile(hierarchyDb)

        // parse the Prints Output into Match instances
        File printsFile = new File(printsOutput)
        Map<String, Map<String, Match>> hits = new LinkedHashMap<>()  // <protein ID <Model Acc, Match>>
        Map<String, String> name2accession = new LinkedHashMap<>()  // <motifName, motifAccession>
        String queryAccession = null  // protein seq ID
        Map<String, Match> thisProteinsMatchesMap = new LinkedHashMap<>()  // modelName: Match()
        Map<String, Set<Match>> allThisProteinMatches = new LinkedHashMap<>()  // all matches for this protein

        printsFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Sn;")) { // Start the new protein: Get the query sequence id
                    queryAccession = line.replaceFirst("Sn; ", "").trim()
                    name2accession = new LinkedHashMap<>()
                    thisProteinsMatchesMap = new LinkedHashMap<>()
                }
                else if (line.startsWith("1TBH")) {
                    // Line: 1TBH 4DISULPHCORE    1.4e-07        4-disulphide core signature       PR00003
                    def matcher1TBH = line =~ ~/^1TBH\s+(\S+)\s+(\S+)\s+.*?\s+(PR\d+)\s+$/
                    if (matcher1TBH.find()) {
                        String modelName = matcher1TBH.group(1).trim()
                        final Double evalue = matcher1TBH.group(2).trim() as Double
                        String modelId = matcher1TBH.group(3).trim()
                        Match match = new Match(modelName)
                        match.evalue = evalue
                        if (!thisProteinsMatchesMap.containsKey(modelName)) {thisProteinsMatchesMap.put(modelName, match)}
                        name2accession.put(modelName, modelId)
                    }
                }
                else if (line.startsWithAny("2TBH", "2TBN")) {  // Retrieve the graphScan value
                    // Line: 2TBH|N  modelId  NumMotifs  SumId  AveId  ProfScore  Ppvalue  Evalue  GraphScan
                    def matcher2TB = line =~ ~/^2TB[H|N]\s+(\S+)\s+\d+\s+of\s+\d+\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\s+.*?\s+.*?\s+([\.iI]+)\s+$/
                    if (matcher2TB.find()) {
                        String modelName = matcher2TB.group(1).trim()
                        final String graphScan = matcher2TB.group(2).trim()
                        if (thisProteinsMatchesMap.containsKey(modelName) && hierarchyMap.containsKey(modelName)) {  // Has to be in the hierarchyDB
                            thisProteinsMatchesMap[modelName].graphScan = graphScan
                        }
                    }
                }
                else if (line.startsWithAny("3TBH", "3TBN")) {  // Used to identify match locations
                    // For the post processing create one Match per location, so each Match obj has one Location obj
                    // Line: 3TBH|N modelId NoOfMotifs IdScore PfScore Pvalue Sequence Len Low Pos High
                    def matcher3TB = line =~ ~/^3TB[H|N]\s+(.+?)\s+(\d+)\s+of\s+(\d+)\s+([\d\.]+)\s+\d+\s+(.+?)\s+([\w#]+)\s+(\d+)\s+\d+\s+(-?\d+)\s*\d\s*$/
                    if (matcher3TB.find()) {
                        String modelName = matcher3TB.group(1).trim()
                        int motifNumber = matcher3TB.group(2).trim() as int  // number of this motif 'X' of Y
                        final Double score = matcher3TB.group(4).trim() as Double
                        final Double pvalue = matcher3TB.group(5).trim() as Double
                        final String motifSequence = matcher3TB.group(6).trim()
                        final int motifLength = matcher3TB.group(7).trim() as int

                        int locationStart = matcher3TB.group(8).trim() as int
                        if (locationStart < 1) {locationStart = 1}

                        int locationEnd = locationStart + motifLength - 1
                        if (motifSequence.endsWith("#")) { // it overhangs the protein seq so adjust locationEnd
                            int motifSeqLength = motifSequence.length()
                            int indexCheck = motifSeqLength - 1
                            while (motifSequence[indexCheck] == "#") {indexCheck -= 1}
                            locationEnd = locationEnd - (motifSeqLength - indexCheck) + 1
                        }

                        Match newMatch = thisProteinsMatchesMap[modelName]
                        Location location = new Location(locationStart, locationEnd, pvalue, score, motifNumber)
                        newMatch.addLocation(location)

                        if (!allThisProteinMatches.containsKey(modelName)) {
                            allThisProteinMatches.put(modelName, new HashSet<Match>())
                        }
                        allThisProteinMatches[modelName].add(newMatch)
                    }

                }
                else if (line.startsWith("3TBF") && !allThisProteinMatches.isEmpty()) { // parse the matches for this protein
                    hits = filterProteinMatches(queryAccession, thisProteinsMatches, hierarchyMap, name2accession, hits)
                }
            }
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
            if (selectMatch(match, modelId, hierarchyEntry, HierarchtEntryIdLimitation)) {
                if (!hits.containsKey(queryAccession)) {hits.put(queryAccession, new LinkedHashMap<>())}
                hits[queryAccession].put(match.modelAccession, match)
                // passed the filter and may have its own hierarchy so update the hierarchy limitation
                if (hierarchyEntry.siblingsIds.length < HierarchtEntryIdLimitation.size()) {
                    HierarchtEntryIdLimitation = hierarchyEntry.siblingsIds
                }
            }
        }

        /// Continue implenting: https://github.com/ebi-pf-team/interproscan/blob/13c095e3c6b6784c98f9769c6a07d122cfa160d1/core/business/src/main/java/uk/ac/ebi/interpro/scan/business/postprocessing/prints/PrintsPostProcessing.java#L226

        return hits
    }

    static sortMatches(List<Match> matches) { // This comparator is CRITICAL to the working of PRINTS post-processing
        return matches.sort { matchA, matchB ->
            int evalueComparison = matchA.evalue <=> matchB.evalue
            if (evalueComparison != 0) return evalueComparison

            int modelAccessionComparison = matchA.modelAccession <=> matchB.modelAccession
            if (modelAccessionComparison != 0) return modelAccessionComparison

            int motifNumberComparison = matchA.locations[0].motifNumber <=> matchB.locations[0].motifNumber
            if (motifNumberComparison != 0) return motifNumberComparison

            int startLocationComparison = matchA.locations[0].start <=> matchB.locations[0].start
            if (startLocationComparison != 0) return startLocationComparison

            return matchA.locations[0].end <=> matchB.locations[0].end
        }
    }

    static selectMatch(Match match, String motifId, HierarchyEntry hierarchy, Set<String> hierarchyEntryIdLimitation) {
        // if a domain: PASS - there is no hierarchy to deal with
        if (hierarchy.isDomain) {return true}
        // if the previous model limited the filtering to its siblings. If the first model, this is all models in HierarchyDB
        if (hierarchyEntryIdLimitation.contains(motifId)) {return true}
        return false
    }
}

class HierarchyEntry implements Serializable { // represents a row in the HierarchyDB
    String modelId
    String modelAccession
    Double evalueCutoff
    int minMotifCount
    boolean isDomain = false
    String[] siblingsIds = []
    
    HierarchyEntry(String modelId, String modelAccession, Double evalueCutoff, int minMotifCount) {
        this.modelId = modelId
        this.modelAccession = modelAccession
        this.evalueCutoff = evalueCutoff
        this.minMotifCount = minMotifCount
    }

    static parseHierarchyDbFile(String hierarchyDb) {  // parser the hierarchyDb in to a map of ModelID: hierarchy entry
        Map<String, HierarchyEntry> hierarchyMap = new LinkedHashMap<>()
        File hierarchyFile = new File(hierarchyDb)
        if (!hierarchyFile.exists()){
            System.out.println("Could not find Hierarchy DB for PRINTS")
            System.exit 1
        }
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
        return hierarchyMap
    }

    void changeDomainStatus(boolean isDomain) {
        this.isDomain = isDomain
    }

    void addSiblings(String[] siblingIds) {
        this.siblingsIds = siblingsIds
    }
}
