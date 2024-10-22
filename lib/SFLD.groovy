class SFLD {
    static parseOutput (
            String sfldPostProcessOutput,
            String sfldHierarchyDb
    ) {
        // get the hierarchy data Map <modelAcc, HierarchyInformation obj>
        Map<String, HierarchyInformation> hierarchyInformation = parseHierarchyFile(sfldHierarchyDb)

        Map<String, Map<String, Match>> hits = new LinkedHashMap<>()  // <protein Id, <Model, Match>>

        File sfldPpOut = new File(sfldPostProcessOutput)
        sfldPpOut.withReader{ reader ->
            String line
            String proteinAccession
            String modelAccession
            String siteResidues
            String siteDescription
            Map<String, Set<Match>> thisProteinsMatches = new LinkedHashMap<>()  // <modelAcc, Set<Matches>>
            Map<String, Set<Site>> thisProteinsSites = new LinkedHashMap<>()     // <modelAcc, Set<Sites>>

            while ((line = reader.readLine()) != null) {
                def seqIdMatcher = line =~ ~/^Sequence:\s+(\S+)$/
                if (seqIdMatcher.find()) {
                    // Process the last protein
                    if (!thisProteinsMatches.isEmpty()) {
                        Map<String, Set<Match>> filteredProteinMatches = filterProteinHits(
                                thisProteinsMatches, hierarchyInformation
                        )
                        Map<String, Set<Match>> filteredProteinMatchesWithSites = filterAndAddProteinSites(
                                filteredProteinMatches, thisProteinsSites
                        )
                        hits.put(proteinAccession, filteredProteinMatchesWithSites)
                    }
                    
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
            Map<String, Set<Match>> filteredProteinMatches = filterProteinHits(
                    thisProteinsMatches, hierarchyInformation
            )
            Map<String, Set<Match>> filteredProteinMatchesWithSites = filterAndAddProteinSites(
                    filteredProteinMatches, thisProteinsSites
            )
            hits.put(proteinAccession, filteredProteinMatchesWithSites)

        }

        return hits
    }

    static parseHierarchyFile(String hierarcyDbPath) { // Map <modelAcc: HierarchyInformation>
        Map<String, HierarchyInformation> hierarchyInformation = new LinkedHashMap<>()
        File hierachyDbFile = new File(hierarcyDbPath)
        if (!hierachyDbFile.exists()) {
            System.out.println("Could not find hierarchy file for SFLD")
            System.exit 1
        }
        hierachyDbFile.withReader{ reader ->
            String line
            while ((line = reader.readLine()) != null) {
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

    static filterProteinHits (Map<String, Set<Match>> rawMatces, Map<String, HierarchyInformation> hierarchyInformationMap) {
        // 1. Remove domains that overlap
        Map<String, Set<Match>> firstFilteredMatches = resolveOverlappingMatches(rawMatces, hierarchyInformationMap)
        // 2. Add matches for parents defined in the hierarchy db
        Map<String, Set<Match>> secondFilteredMatches = addParentMatches(firstFilteredMatches, hierarchyInformationMap)
        // 3. Remove duplicated hits
        Map<String, Set<Match>> thirdFilteredMatches = resolveDuplicateMatches(secondFilteredMatches)
        return thirdFilteredMatches
    }

    static resolveOverlappingMatches (
            Map<String, Set<Match>> currentMatches,
            Map<String, HierarchyInformation> hierarchyInformation
    ) {
        Set<Match> allRawMatches = currentMatches.values().collect().flatten().toSet()
        Map<String, Set<Match>> nonOverlappingMatches = new LinkedHashMap<>()

        for (String modelAccession: currentMatches.keySet()) {
            if (currentMatches[modelAccession] == null) {
                continue
            }
            if (modelAccession.contains("SFLDF") || currentMatches[modelAccession].size() == 1) {
                if (!nonOverlappingMatches.containsKey(modelAccession)) {
                    nonOverlappingMatches.put(modelAccession, [] as Set<Match>)
                }
                nonOverlappingMatches[modelAccession].addAll(currentMatches[modelAccession])
                continue
            }

            for (Match match: currentMatches[modelAccession]) {
                boolean overlaps = false
                for (Match otherMatch: allRawMatches) {
                    String otherMatchParents = hierarchyInformation.get(otherMatch.modelAccession)
                    if (match == otherMatch || modelAccession == otherMatch.modelAccession) {
                        continue
                    }
                    if ((match.locations[0].start < otherMatch.locations[0].end ||
                                    otherMatch.locations[0].start < match.locations[0].end)
                                    && otherMatchParents.contains(modelAccession)) {
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
        Map<String, Set<Match>> filteredMatches = currentMatches  // we will add the parent matches to these
        Set<Match> seqMatches = [] as Set<Match>

        for (String modelAccession: currentMatches) {
            HierarchyInformation hierarchyRecord = hierarchyInformationMap.get(modelAccession)
            if (hierarchyRecord != null ) {
                if (hierarchyRecord.parents != null && hierarchyRecord.parents.size() > 0) {

                    for (Match match: currentMatches[modelAccession]) {
                        Set<Match> parentMatches = getParents(match, hierarchyRecord.parents)
                        boolean parentContainsSeqMatch = false
                        boolean seqMatchContainsParentMatch = false
                        Match matchToRemove = null
                        for (Match parentMatch: parentMatches) {
                            for (Match seqMatch: seqMatches) {
                                if ( parentMatch.locations[0].start <= seqMatch.locations[0].start &&
                                    parentMatch.locations[0].end >= seqMatch.locations[0].end) {
                                    parentContainsSeqMatch = true
                                    matchToRemove = parentMatch
                                }
                                if ( seqMatch.locations[0].start <= parentMatch.locations[0].start &&
                                    seqMatch.locations[0].end >= parentMatch.locations[0].end) {
                                    seqMatchContainsParentMatch = true
                                    break
                                }
                            }
                            if (parentContainsSeqMatch) {seqMatches.remove(matchToRemove)}
                            if (seqMatchContainsParentMatch) {seqMatches.add(parentMatch)}
                        } // end of parent Match
                    } // end of Match
                }
            }
        } // end of modelAccession

        for (Match match: seqMatches) {
            if (!filteredMatches.containsKey(match.modelAccession)) {
                filteredMatches.put(match.modelAccession, [] as Set<Match>)
            }
            filteredMatches[match.modelAccession].add(match)
        }

        return filteredMatches
    }

    static getParents(Match childMatch, Set<String> parents) {
        // Called in addParentMatches
        // Build parent matches from the child match
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

    static resolveDuplicateMatches(
            Map<String, Set<Match>> currentMatches
    ) {
        Map<String, Set<Match>> duplicateFreeMatches = new LinkedHashMap<>()
        Set<Match> allMatches = currentMatches.values().collect().flatten().toSet()

        for (String modelAccession: currentMatches.keySet()) {
            boolean duplicate = false
            for (Match match: currentMatches[modelAccession]) {
                for (Match otherMatch: allMatches) {
                    if (match == otherMatch || !otherMatch.modelAccession == modelAccession) { continue }
                    if (match.locations[0].start == otherMatch.locations[0].start &&
                        match.locations[0].end == otherMatch.locations[0].end &&
                        match.locations[0].score <= otherMatch.locations[0].score) {
                        duplicate = true
                        break
                    }
                } //end of Other match
                if (duplicate) {
                    allMatches.remove(match)
                }
                else {
                    if (!duplicateFreeMatches.containsKey(modelAccession)) {
                        duplicateFreeMatches.put(modelAccession, [] as Set<Match>)
                    }
                    duplicateFreeMatches[modelAccession].add(match)
                }
            } // end for match
        } // end for model Accession
        return duplicateFreeMatches
    }

    static filterAndAddProteinSites (Map<String, Set<Match>> filteredMatches, Map<String, Set<Site>> rawSites) {
        // add the sites retrieved from SFLD binary output file to the corresponding filtered match (if one exists)
        Map<String, Set<Match>> filteredMatchesWithSites = new LinkedHashMap<>()

        for (String modelAccession: filteredMatches.keySet()) {
            if (!filteredMatchesWithSites.containsKey(modelAccession)) {
                filteredMatchesWithSites.put(modelAccession, [] as Set<Match>)
            }

            for (Match match: filteredMatches[modelAccession]) {
                if (rawSites.containsKey(modelAccession)) {
                    for (Site site: rawSites[modelAccession]) {
                        match.addSite(site) // will automatically add the site to the correct location
                    }
                }
                // gone though all sites (if any), so add the match to the new map
                filteredMatchesWithSites[modelAccession].add(match)
            }
        }
        return filteredMatchesWithSites
    }
}


class HierarchyInformation {
    String childModelAcc
    Set<String> parents = [] as Set<String>

    HierarchyInformation(String childModelAcc) {
        this.childModelAcc = childModelAcc
    }

    void addParent(String parentModelAcc) {
        this.parents.add(parentModelAcc)
    }
}