process RUN_PRINTS {
    label     'mem_medium'
    label     'time_short'
    container 'interpro/fingerprintscan:3.597.ebiftp'

    input:
    tuple val(meta), path(fasta)
    path pval

    output:
    tuple val(meta), path("prints_output")

    script:
    """
    fingerPRINTScan \
        ${pval} \
        ${fasta} \
        -e 0.0001 -d 10 -E 257043 84355444 -fj -o 15 > prints_output
    """
}

process PARSE_PRINTS {
    label    'mem_low'
    label    'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(prints_output)
    val hierarchydb

    output:
    tuple val(meta), path("prints.json")

    exec:
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("PRINTS", null)
    // Build up a map of the Model ID to fingerprint hierarchies
    def hierarchyMap = uk.ac.ebi.interpro.FingerPrint.HierarchyEntry.parseHierarchyDbFile(hierarchydb)

    // Parse the prints output into simple raw prints matches
    // Each location is represented by its own Print object
    def queryAccession = null                        // protein seq ID
    def thisProteinsMatches = [:]  // modelName: FingerPrint()
    def rawMatches = [:]     // <protein ID <FingerPrint()>>
    prints_output.withReader { reader ->
        reader.eachLine { line ->
            if (line.startsWith("Sn;")) { // Start the new protein: Get the query sequence id
                queryAccession = line.split(/\s+/)[1]
                thisProteinsMatches = [:]
            } else if (line.startsWith("1TBH")) {
                // Line: 1TBH 4DISULPHCORE    1.4e-07        4-disulphide core signature       PR00003
                // we retrieve the model description during the XREFs stage of IPS6
                def lineData1TBH = line.split(/\s+/)
                assert lineData1TBH.length >= 5

                def modelName = lineData1TBH[1]
                assert hierarchyMap[modelName] != null  // the model must be in the Hierarchy DB
                def evalue = lineData1TBH[2] as double
                def modelId = lineData1TBH[-1]

                def printMatch = new uk.ac.ebi.interpro.FingerPrint(modelName, modelId, evalue)
                thisProteinsMatches.computeIfAbsent(modelName, {printMatch})
            } else if (line.startsWithAny("2TBH", "2TBN")) {
                // Line: 2TBH|N  modelId  NumMotifs  SumId  AveId  ProfScore  Ppvalue  Evalue  graphscan
                // Retrieve the graphscan value
                def lineData2TBHN = line.split(/\s+/)
                assert lineData2TBHN.length == 11
                def modelName = lineData2TBHN[1]
                def graphscan = lineData2TBHN[-1]
                if (thisProteinsMatches.containsKey(modelName)) {
                    thisProteinsMatches[modelName].graphscan = graphscan
                }
            } else if (line.startsWithAny("3TBH", "3TBN")) {
                // For the post processing create one Match per location, so each Match obj has one Location obj
                // Line: 3TBH|N modelId NoOfMotifs IdScore PfScore Pvalue Sequence Len Low Pos High
                def lineData3TBHN = line.split(/\s+/)
                assert lineData3TBHN.length == 13

                def modelName = lineData3TBHN[1]
                if (thisProteinsMatches.containsKey(modelName)) {
                    def motifNumber = lineData3TBHN[2] as int  // number of this motif 'X' of Y -- take the X
                    def motifCount = lineData3TBHN[4] as int // X of 'Y' -- take the Y
                    def score = lineData3TBHN[5] as double
                    def pvalue = lineData3TBHN[7] as double
                    def motifSequence = lineData3TBHN[8]
                    def motifLength = lineData3TBHN[9] as int

                    def locationStartString
                    if (lineData3TBHN[11].length() > 5) {
                        // A starting position that is more than 5 figures merges into the high column
                        locationStartString = lineData3TBHN[11].take(5)
                    } else {
                        locationStartString = lineData3TBHN[11]
                    }

                    def locationStart = locationStartString as int
                    def locationEnd = locationStart + motifLength - 1

                    // Adjust the locationStart if the motif starts before the sequence
                    locationStart = Math.max(1, locationStart) 

                    if (motifSequence.endsWith("#")) { 
                        // it overhangs the protein seq so adjust locationEnd
                        def trailingHashCount = motifSequence.length() - motifSequence.replaceFirst(/#+$/, "").length()
                        locationEnd -= trailingHashCount
                    }

                    def matchLocation = uk.ac.ebi.interpro.FingerPrint.buildLocationMatch(
                            thisProteinsMatches[modelName],
                            locationStart,
                            locationEnd,
                            pvalue,
                            score,
                            motifNumber,
                            motifCount
                    )
                    rawMatches
                        .computeIfAbsent(queryAccession, {[]})
                        .add(matchLocation)
                }
            }        
        }
    }

    /* Filter the prints matches
    General overview:
    Matches to "Domain" PRINTS models can match outside of the hierarchy constraints, but apply the
    hierarchy constraints to all other matches.
    Filter matches:
    1. evalue <= cutOff defined in hierarchyDb
    2. motif >= min number of motifs defined in the hierarchy Db
    3. order by evalue (best first)
    4. If a domain pass
    5. Apply the hierarchy constraints to see if the match passes
    */
    def matches = [:]
    rawMatches.each { proteinAccession, proteinMatches ->
        def sortedMatches = sortMatches(proteinMatches)
        def filteredMatches = []
        def currentModelAcc = null
        def motifMatchesForCurrentModel = []
        def currentMatchesPass = true
        def passed = false
        def currentHierarchyEntry = null
        def hierarchyEntryIdLimitation = hierarchyMap.keySet() // initialise with all models

        sortedMatches.each { rawMatch ->
            if (currentModelAcc == null || currentModelAcc != rawMatch.modelName) {
                // just started or moved onto a match for a different model
                if (currentModelAcc != null && currentMatchesPass) {
                    passed = selectMatches(
                        motifMatchesForCurrentModel, currentModelAcc,
                        currentHierarchyEntry, hierarchyEntryIdLimitation
                    )
                    if (passed) {
                        filteredMatches.addAll(motifMatchesForCurrentModel)
                        if (currentHierarchyEntry.siblingsIds.size() < hierarchyEntryIdLimitation.size() 
                                && !currentHierarchyEntry.isDomain) {
                            hierarchyEntryIdLimitation = currentHierarchyEntry.siblingsIds.toSet()
                        }
                    }
                }

                // reset the values
                currentMatchesPass = true
                motifMatchesForCurrentModel = []
                currentModelAcc = rawMatch.modelName
                assert hierarchyMap[currentModelAcc] != null
                currentHierarchyEntry = hierarchyMap[currentModelAcc]
            }
            if (currentMatchesPass) { currentMatchesPass = rawMatch.evalue <= currentHierarchyEntry.evalueCutoff }
            if (currentMatchesPass) { motifMatchesForCurrentModel.add(rawMatch) }
        }
        // parse the last matches
        passed = selectMatches(
             motifMatchesForCurrentModel, currentModelAcc,
             currentHierarchyEntry, hierarchyEntryIdLimitation
        )
        if (passed) {
            filteredMatches.addAll(motifMatchesForCurrentModel)
        }

        // add the filteredMatches to matches
        if (!filteredMatches.isEmpty()) {
            def finalMatches = matches.computeIfAbsent(proteinAccession, { [:] })
            filteredMatches.each { filteredMatch ->
                // the modelName has been used up to this point, but we need to convert to the model ID
                def match = finalMatches.computeIfAbsent(
                    filteredMatch.modelId,
                    {
                        new uk.ac.ebi.interpro.Match(
                            filteredMatch.modelId, 
                            filteredMatch.evalue, 
                            filteredMatch.graphscan, 
                            new uk.ac.ebi.interpro.Signature(filteredMatch.modelId, library))
                    }
                )
                def location = new uk.ac.ebi.interpro.Location(
                    filteredMatch.locationStart,
                    filteredMatch.locationEnd,
                    filteredMatch.pvalue,
                    filteredMatch.score,
                    filteredMatch.motifNumber
                )
                match.addLocation(location)
            }
        }
    }

    def filepath = task.workDir.resolve("prints.json")
    filepath.text  = groovy.json.JsonOutput.toJson(matches)
}

def sortMatches(matches) {
    // This comparator is CRITICAL to the working of PRINTS post-processing
    return matches.sort { matchA, matchB ->
        def evalueComparison = matchA.evalue <=> matchB.evalue
        if (evalueComparison != 0) return evalueComparison

        def modelAccessionComparison = matchA.modelId <=> matchB.modelId
        if (modelAccessionComparison != 0) return modelAccessionComparison

        def motifNumberComparison = matchA.motifNumber <=> matchB.motifNumber
        if (motifNumberComparison != 0) return motifNumberComparison

        def startLocationComparison = matchA.locationStart <=> matchB.locationStart
        if (startLocationComparison != 0) return startLocationComparison

        return matchA.locationEnd <=> matchB.locationEnd
    }
}

def selectMatches(motifMatchesForCurrentModel, modelName, hierarchy, hierarchyEntryIdLimitation) {
    // Belt and braces check
    if (motifMatchesForCurrentModel.size() == 0) { return false }
    // Check that enough motifs for the current model passed the previous filtering criteria
    if (motifMatchesForCurrentModel.size() < hierarchy.minMotifCount) { return false }
    // If the model is a domain: PASS - there is no hierarchy to deal with
    if (hierarchy.isDomain) { return true }
    // Check the model meets the hierarchy filtering criteria
    if (!hierarchyEntryIdLimitation.contains(modelName)) { return false }
    return true
}
