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
        Set<Match> allThisProteinMatches = [] as Set<Match>  // all matches for this protein
        Set<Match> proteinMatches = [] as Set<Match>

        printsFile.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                if (line.startsWith("Sn;")) { // Start the new protein: Get the query sequence id
                    queryAccession = line.replaceFirst("Sn; ", "").split(" ")[0].trim()
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

                        if (thisProteinsMatchesMap.containsKey(modelName)) {
                            int motifNumber = matcher3TB.group(2).trim() as int  // number of this motif 'X' of Y
                            final Double score = matcher3TB.group(4).trim() as Double
                            final Double pvalue = matcher3TB.group(5).trim() as Double
                            final String motifSequence = matcher3TB.group(6).trim()
                            final int motifLength = matcher3TB.group(7).trim() as int

                            int locationStart = matcher3TB.group(8).trim() as int
                            if (locationStart < 1) {
                                locationStart = 1
                            }

                            int locationEnd = locationStart + motifLength - 1
                            if (motifSequence.endsWith("#")) { // it overhangs the protein seq so adjust locationEnd
                                int motifSeqLength = motifSequence.length()
                                int indexCheck = motifSeqLength - 1
                                while (motifSequence[indexCheck] == "#") {
                                    indexCheck -= 1
                                }
                                locationEnd = locationEnd - (motifSeqLength - indexCheck) + 1
                            }

                            Match newMatch = thisProteinsMatchesMap[modelName]
                            Location location = new Location(locationStart, locationEnd, pvalue, score, motifNumber)
                            newMatch.addLocation(location)
                            allThisProteinMatches.add(newMatch)
                        }
                    }

                }
                else if (line.startsWith("3TBF") && !allThisProteinMatches.isEmpty()) { // parse the matches for this protein
                    System.out.println("${queryAccession} Matches: ${allThisProteinMatches.size()}")
                    proteinMatches = filterProteinMatches(queryAccession, allThisProteinMatches, hierarchyMap, name2accession, hits)
                    System.out.println("${queryAccession} Filtered Matches: ${proteinMatches.size()}")
                    hits = storeMatches(hits, proteinMatches as List<Match>, queryAccession, name2accession)
                    name2accession = new LinkedHashMap<>()  // <motifName, motifAccession>
                    queryAccession = null  // protein seq ID
                    thisProteinsMatchesMap = new LinkedHashMap<>()  // modelName: Match()
                    allThisProteinMatches = [] as Set<Match>  // all matches for this protein
                    proteinMatches = [] as Set<Match>
                }
            }
        }
        return hits
    }

    static filterProteinMatches(
            String queryAccession,
            Set<Match> thisProteinsMatches,
            Map<String, HierarchyEntry> hierarchyMap,
            Map<String, String> name2accession,
            Map<String, Map<String, Match>> hits
    ) {
        /* General overview:
        Matches to "Domain" PRINTS models can match outside of the hierarchy constraints, but apply the
        hierarchy constraints to all other matches.
        Filter matches:
        1. evalue <= cutOff defined in hierarchyDb
        2. motif >= min number of motifs defined in the hierarchy Db
        3. order by evalue (best first)
        4. If a domain pass
        5. Apply the hierarchy constraints to see if the match passes
         */
        List<Match> sortedMatches = sortMatches(thisProteinsMatches as List<Match>)
        System.out.println("Inside filterProteinMatches: ${sortedMatches.size()} sortedMatches")
        List<Match> filteredMatches = [] as List<Match>

        String currentModelAcc = null
        List<Match> motifMatchesForCurrentModel = [] as List<Match>
        boolean currentMatchesPass = true
        boolean passed = false
        HierarchyEntry currentHierarchyEntry = null
        Set<String> hierarchtEntryIdLimitation = hierarchyMap.keySet() // initialise with all model

        for (Match match in sortedMatches) {
            if (currentModelAcc == null || currentModelAcc != match.modelAccession) {
                // just started or moved onto the matches for the next model
                if (currentModelAcc != null && currentMatchesPass && motifMatchesForCurrentModel.size() > 0) {
                    passed = selectMatches(
                            motifMatchesForCurrentModel, currentModelAcc,
                            currentHierarchyEntry, hierarchtEntryIdLimitation
                    )
                    if (passed) {
                        System.out.println("Passes, ${motifMatchesForCurrentModel}")
                        filteredMatches.addAll(motifMatchesForCurrentModel)
                        if (currentHierarchyEntry.siblingsIds.length < hierarchtEntryIdLimitation.size()) {
                            hierarchtEntryIdLimitation = currentHierarchyEntry.siblingsIds
                        }
                    }
                }
                // reset values
                currentMatchesPass = true
                motifMatchesForCurrentModel.clear()
                currentModelAcc = match.modelAccession
                currentHierarchyEntry = hierarchyMap[match.modelAccession]
            }

            if (currentMatchesPass) {currentMatchesPass = match.evalue <= currentHierarchyEntry.evalueCutoff}
            if (currentMatchesPass) {motifMatchesForCurrentModel.add(match)}
        }

        // parse the last matches
        if (motifMatchesForCurrentModel.size() > 0) {
            if (selectMatches(
                    motifMatchesForCurrentModel, currentModelAcc,
                    currentHierarchyEntry, hierarchtEntryIdLimitation
            )) {filteredMatches.addAll(motifMatchesForCurrentModel)}
        }

        return filteredMatches
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

    static selectMatches(
        List<Match> motifMatchesForCurrentModel,
        String modelName,
        HierarchyEntry hierarchy,
        Set<String> hierarchyEntryIdLimitation
    ) {
        System.out.println("len hier ${hierarchyEntryIdLimitation.size()}")
        System.out.println("modelName ${modelName} -- ${hierarchyEntryIdLimitation[0]}")

        // Check that enough motifs for the current model passed the previous filtering criteria
        if (motifMatchesForCurrentModel.size() < hierarchy.minMotifCount) {
            System.out.println("*XX* too few motifs")
            return false }
        // If the model is a domain: PASS - there is no hierarchy to deal with
        if (hierarchy.isDomain) {
            System.out.println("** isDomain")
            return true }
        // Check the model meets the hierarchy filtering criteria
        if (hierarchyEntryIdLimitation.contains(modelName)) {
            System.out.println("** inHierarchy")
            return true }
        System.out.println("*XX* notInHierarchy")
        return false
    }

    static storeMatches(
            Map<String, Map<String, Match>> hits,
            List<Match> proteinMatches,
            String proteinId,
            Map<String, String> name2accession
    ) {
        if ( !hits.containsKey(proteinId) ) { hits.put(proteinId, new LinkedHashMap<String, Match>()) }

        for (Match match: proteinMatches) {
            // the modelName has been used up to this point, but we need to convert to the correct model Accession/Id
            String trueModelAccession = name2accession.get(match.modelAccession)
            if (!hits[proteinId].containsKey(trueModelAccession)) {
                Match newMatch = new Match(trueModelAccession, match.evalue, match.graphScan)
                hits[proteinId].put(trueModelAccession, newMatch)
            }
            Location newLocation = new Location(
                    match.locations[0].start,
                    match.locations[0].end,
                    match.locations[0].pvalue,
                    match.locations[0].score,
                    match.locations[0].motifNumber,
            )
            hits[proteinId][trueModelAccession].addLocation(newLocation)
        }
        return hits
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
