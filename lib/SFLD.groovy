class SFLD {
    static parseOutput (
            String sfldPostProcessOutput,
            String sfldHierarchyDb
    ) {
        // get the hierarchy data Map <modelAcc, HierarchyInformation obj>
        Map<String, HierarchyInformation> hierarchyInformation = HierarchyInformation.parseHierarchyFile(sfldHierarchyDb)

        Map<String, Map<String, Match>> hits = new LinkedHashMap<>()  // <ModelAcc <protein ID, Match>>

        File sfldPpOut = new File(sfldPostProcessOutput)
        sfldPpOut.withReader{ reader ->
            String line
            String proteinAccession
            String modelAccession
            String siteResidues
            String siteDescription
            Map<String, Set<Match>> thisProteinsMatches = new LinkedHashMap<>()  // <modelAcc, Set<Matches>>
            Map<String, Set<Site>> thisProteinsSites = new LinkedHashMap<>()    // <modelAcc, Set<Sites>>

            while ((line = reader.readLine() != null)) {
                def seqIdMatcher = line =~ ~/^Sequence:\s+(.+)$/
                if (seqIdMatcher.find()) {
                    // Process the last protein
                    //
                    //
                    
                    // Initialise the new protein
                    proteinAccession = seqIdMatcher.group(1).trim()
                    thisProteinsMatches = new LinkedHashMap<>()  // <modelAcc, Set<Matches>>
                    thisProteinsSites = new LinkedHashMap<>()    // <modelAcc, Set<Sites>>
                    continue
                }
                
                // Retrieve domain data
                // e.g. SFLDF00315	0.000e+00	6.097e+02	0.300	1	443	609.600	1	447	1	448	0.000e+00	0.000e+00	0.990	0.300
                // modelAcc seqEvalue seqScore seqBias hmmFrom hmmTo domScore aliFrom aliTo envFrom envTo domCevalue DomIevalue acc domBias
                def domainMatcher = line =~ ~/^(SFLD\D+\d+)\s(\d+\.\d+e[+\-]\d+)\s+(\d+\.\d+e[+-]\d+)\s+\d\.?\d*?\s+(\d+)\s+(\d+)\s+(\d+\.?\d*?)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+\d+\.\d+e[+-]\d+\s+(\d+\.\d+e[+-]\d+)\s+\d+\.\d+\s(\d+\.\d+)\s?$/
                if (domainMatcher.find()) {
                    // Define the Match instance
                    modelAccession = domainMatcher.group(1).trim()
                    Double seqEvalue = domainMatcher.group(2).trim() as Double
                    Double seqScore = domainMatcher.group(3).trim() as Double
                    Match match = new Match(modelAccession, seqEvalue, seqScore)
                    
                    // Define the location for this hit
                    int start = domainMatcher.group(7).trim() as int    // aliFrom
                    int end = domainMatcher.group(8).trim() as int      // aliTo
                    Integer hmmStart = domainMatcher.group(4).trim() as Integer
                    Integer hmmEnd = domainMatcher.group(5).trim() as Integer
                    Double domainEvalue = domainMatcher.group(11) as Double
                    Double domainScore = domainMatcher.group(6).trim() as Double
                    Integer envStart = domainMatcher.group(9).trim() as Integer
                    Integer envEnd = domainMatcher.group(10).trim() as Integer
                    Double domainBias = domainMatcher.group(12).trim() as Double
                    List<LocationFragment> fragment = [new  LocationFragment(start, end, "CONTINUOUS")]
                    // hmmLength is retrieved from the hmmersearch.dtbl file and added to the match later
                    Location location = new Location(
                            start, end, hmmStart, hmmEnd, envStart, envEnd, domainEvalue, domainScore, domainBias, fragment
                    )
                    match.addLocation(location)
                    
                    if (!thisProteinsMatches.containsKey(modelAccession)) {
                        thisProteinsMatches.put(modelAccession, [] as Set<Match>)
                    }
                    thisProteinsMatches[modelAccession].add(match)
                    
                    continue
                }
                
                // Retrieve Site data
                // e.g. SFLDG01386 C129,C133,C136 Binds [4Fe-4S]-AdoMet cluster
                def siteMatcher = line =~ ~/^(SFLD\D+\d+)\s+((?:\D\d+,?)+)\s+(.+)$/
                if (siteMatcher.find()) {
                    modelAccession = siteMatcher.group(1).trim()
                    siteResidues = siteMatcher.group(2).trim()
                    siteDescription = siteMatcher.group(3).trim()
                    Site site = new Site(siteDescription, siteResidues)
                    if (!thisProteinsSites.containsKey(modelAccession)) {
                        thisProteinsSites.put(modelAccession, [] as Set<Site>)
                    }
                    thisProteinsSites[modelAccession].add(site)
                }
            }

            // Parse the last protein
            thisProteinsMatches = resolveOverlappingMatches(thisProteinsMatches, hierarchyInformation)
            thisProteinsMatches = filterMatches(thisProteinsMatches, hierarchyInformation)

        }

        return hits
    }

    static filterProteinHits (Map<String, Set<Match>> rawMatces, Map<String, HierarchyInformation> hierarchyInformationMap) {
        // 1. Remove domains that overlap
        Map<String, Set<Match>> firstFilteredMatches = resolveOverlappingMatches(rawMatces, hierarchyInformationMap)
        // 2. Add matches for parents defined in the hierarchy db
        Map<String, Set<Match>> secondFilteredMatches = addParentMatches(firstFilteredMatches, hierarchyInformationMap)
        // 3.
    }

    static filterProteinSites (Map<String, Set<Match>> filteredMatches, Map<String, Set<Site>> rawSites) {
        // add the sites retrieved from SFLD binary output file to the corresponding filtered match (if one exists)
    }

    static resolveOverlappingMatches (
            Map<String, Set<Match>> currentMatches,
            Map<String, HierarchyInformation> hierarchyInformation
    ) {
        Set<Match> allRawMatches = currentMatches.values().collect().flatten().toSet()
        Map<String, Set<Match>> nonOverlappingMatches = new LinkedHashMap<>()

        for (String modelAccession: currentMatches) {
            if (modelAccession.contains("SFLDF") || currentMatches[modelAccession].size() == 1) {
                nonOverlappingMatches.addAll(currentMatches[modelAccession])
                continue
            }

            for (Match match: currentMatches[modelAccession]) {
                boolean overlaps = false
                for (Match otherMatch: allRawMatches) {
                    String otherMatchParents = hierarchyInformation.get(otherMatch.modelAccession)
                    if (match == otherMatch || modelAccession == otherMatch.modelAccession) {
                        continue
                    }
                    if (matchesOverlap(match, otherMatch) && otherMatchParents.contains(modelAccession)) {
                        overlaps = true
                        break
                    }
                }
                if (overlaps) { allRawMatches.remove(match) }
                else {
                    if (!nonOverlappingMatches.containsKey(modelAccession)) {
                        nonOverlappingMatches.put(modelAccession, [] as Set<Match>)
                    }
                    nonOverlappingMatches[modelAccession].add(match)
                }
            }
        }

        return nonOverlappingMatches
    }

    static addParentMatches (
            Map<String, Set<Match>> currentMatches,
            Map<String, HierarchyInformation> hierarchyInformationMap
    ) {
        Set<Match> allRawMatches = currentMatches.values().collect().flatten().toSet()
        Map<String, Set<Match>> filteredMatches = currentMatches  // we will add the parent matches to these

        for (String modelAccession: currentMatches) {
            HierarchyInformation hierarchyRecord = hierarchyInformationMap.get(modelAccession)
            if (hierarchyRecord != null ) {
                if (hierarchyRecord.parents != null && hierarchyRecord.parents.size() > 0) {

                    for (Match match: currentMatches[modelAccession]) {
                        parentMatches = getParentMatches(match, hierarchyRecord.parents)

                    }
                }
            }
        }
    }

    static getParentMatches(Match childMatch, Set<String> parents) { // Build parent matches from the child match
        Set<Match> parentMatches = []
        for (String parentModelAcc in parents) {
            if (childMatch.modelAccession != parentModelAcc) {
                Match parentMatch = childMatch.clone()
                parentMatch.modelAccession = parentModelAcc
                parentMatches.add(parentMatch)
            }
        }
        return parentMatches
    }



    static filterMatches (List<Match> currentMatches, HierarchyInformation hierarchyInfo) {
        Set<Match> selectedMatches = []  // matches selected for the current protein sequence
        for (Match match: currentMatches) {
            Set<String> parents = hierarchyInfo.getParents(match.modelAccession, hierarchyInfo)  // Set of model Accs
            Set<Match> parentMatches = []

            if (parents != null && parents.size() > 0) {
                parentMatches = getParentMatches(match, parents)
                boolean parentContainsSelectedMatch = false
                boolean selectedMatchContainsParentMatch = false
                Match matchToRemove = null
                // Compare the currently selected matches to the current parent matches
                for (Match parentMatch: parentMatches) {
                    for (Match selectedMatch: selectedMatches) {
                        // Does the parentMatch match or extend beyond the selectedMatch?
                        if (parentMatch.getLocationStart() <= selectedMatch.getLocationStart() &&
                                parentMatch.getLocationEnd() >= selectedMatch.getLocationEnd()) {
                            parentContainsSelectedMatch = true
                            matchToRemove = selectedMatch
                        }
                        // Does the selectedMatch match or extend beyond the parentMatch?
                        else if (selectedMatch.getLocationStart() <= parentMatch.getLocationStart() &&
                                selectedMatch.getLocationEnd() >= selectedMatch.getLocationEnd()) {
                            selectedMatchContainsParentMatch = true
                        }
                    }
                    // Update the selectedMatches
                    if (parentContainsSelectedMatch) {selectedMatches.remove(matchToRemove)}
                    if (!selectedMatchContainsParentMatch) {selectedMatches.add(parentMatch)}
                }
            }
        }

        // Rewriting this code - https://github.com/ebi-pf-team/interproscan/blob/13c095e3c6b6784c98f9769c6a07d122cfa160d1/core/io/src/main/java/uk/ac/ebi/interpro/scan/io/match/hmmer/hmmer3/SFLDHmmer3MatchParser.java#L433
        if (selectedMatches.size() > 0) {
            Set<Match> allMatches = selectedMatches
            Set<Match> duplicateFreeMatches = []
            for (Match match: selectedMatches) {
                boolean isDuplicate = false
                for (Match otherMatch: selectedMatches) {
                    // if it's itself or for different signatures, SKIP!
                    if (match == otherMatch || match.modelAccession != otherMatch.modelAccession) {
                        continue
                    }
                    // then check the start, end and score
                    if (match.getLocationStart() == otherMatch.getLocationStart() &&
                            match.getLocationEnd() == otherMatch.getLocationEnd() &&
                            match.score <= otherMatch.score) {
                        isDuplicate = true
                        break
                    }
                }
                if (!isDuplicate) {duplicateFreeMatches.add(Match)}
            }
        }
    }

    static matchesOverlap(Match match, Match otherMatch) {
        return match.locations[0].start < otherMatch.locations[0].end || otherMatch.locations[0].start < match.locations[0].end
    }


}


class HierarchyInformation {
    String childModelAcc
    Set<String> parents

    HierarchyInformation(String childModelAcc) {
        this.childModelAcc = childModelAcc
    }

    static addParent(String parentModelAcc) {
        this.parents.add(parentModelAcc)
    }

    static parseHierarchyFile(String hierarcyDbPath) { // Map <modelAcc: HierarchyInformation>
        Map<String, HierarchyInformation> hierarchyInformation = new LinkedHashMap<>()
        File hierachyDbFile = new File(hierarcyDbPath)
        hierachyDbFile.withReader{ reader ->
            String line
            while ((line = reader.readLine() != null)) {
                String[] modelWithParents = line.trim().split(":")
                if (modelWithParents.length >= 2){
                    String modelAccession = modelWithParents[0]
                    HierarchyInformation hierarchyRecord = new HierarchyInformation(modelAccession)
                    for (String parent: modelWithParents[1].split("\\s+")) {
                        if (!parent.trim().isEmpty() || parent.trim() == modelAccession) {
                            hierarchyRecord.addParent(parent.trim())
                        }
                    }
                    hierarchyInformation.put(modelAccession, hierarchyRecord)
                }
            }
        }
        return hierarchyInformation
    }
}