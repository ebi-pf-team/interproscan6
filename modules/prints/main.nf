import groovy.json.JsonOutput
import Prints
import HierarchyEntry

import Match

process RUN_PRINTS {
    label 'medium', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val pval

    output:
    tuple val(meta), path("prints_output")

    script:
    """
    $projectDir/bin/prints/fingerPRINTScan \
        ${dirpath}/${pval} \
        ${fasta} \
        -e 0.0001 -d 10 -E 257043 84355444 -fj -o 15 > prints_output
    """
}

process PARSE_PRINTS {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(prints_output)
    val dirpath
    val hierarchydb

    output:
    tuple val(meta), path("prints.json")

    exec:
    SignatureLibraryRelease library = new SignatureLibraryRelease("PRINTS", null)
    String hierarchyFilePath = "${dirpath.toString()}/${hierarchydb}"
    // Build up a map of the Model ID to fingerprint hierarchies
    Map<String, HierarchyEntry> hierarchyMap = HierarchyEntry.parseHierarchyDbFile(hierarchyFilePath)

    // Parse the prints output into simple raw prints matches
    // Each location is represented by its own Print object
    File printsFile = new File(prints_output.toString())
    String queryAccession = null                   // protein seq ID
    Map<String, Prints> thisProteinsMatches = [:]  // modelName: Prints()
    Map<String, List<Prints>> rawMatches = [:]      // <protein ID <Prints()>>
    printsFile.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null) {

            if (line.startsWith("Sn;")) { // Start the new protein: Get the query sequence id
                queryAccession = line.split(/\s+/)[1]
                thisProteinsMatches = new LinkedHashMap<>()
            }

            else if (line.startsWith("1TBH")) {
                // Line: 1TBH 4DISULPHCORE    1.4e-07        4-disulphide core signature       PR00003
                // we retrieve the model description during the XREFs stage of IPS6
                def lineData1TBH = line.split(/\s+/)
                assert lineData1TBH.length >= 5

                String modelName = lineData1TBH[1]
                assert hierarchyMap[modelName] != null  // the model must be in the Hierarchy DB
                double evalue = lineData1TBH[2] as double
                String modelId = lineData1TBH[-1]

                Prints printMatch = new Prints(modelName, modelId, evalue)
                thisProteinsMatches.computeIfAbsent(modelName, {printMatch})
            }

            else if (line.startsWithAny("2TBH", "2TBN")) {
                // Line: 2TBH|N  modelId  NumMotifs  SumId  AveId  ProfScore  Ppvalue  Evalue  graphscan
                // Retrieve the graphscan value
                def lineData2TBHN = line.split(/\s+/)
                assert lineData2TBHN.length == 11
                String modelName = lineData2TBHN[1]
                String graphscan = lineData2TBHN[-1]
                if (thisProteinsMatches.containsKey(modelName)) {
                    thisProteinsMatches[modelName].graphscan = graphscan
                }
            }

            else if (line.startsWithAny("3TBH", "3TBN")) {
                // For the post processing create one Match per location, so each Match obj has one Location obj
                // Line: 3TBH|N modelId NoOfMotifs IdScore PfScore Pvalue Sequence Len Low Pos High
                def lineData3TBHN = line.split(/\s+/)
                assert lineData3TBHN.length == 13

                String modelName = lineData3TBHN[1]
                if (thisProteinsMatches.containsKey(modelName)) {
                    int motifNumber = lineData3TBHN[2] as int  // number of this motif 'X' of Y -- take the X
                    int motifCount = lineData3TBHN[4] as int // X of 'Y' -- take the Y
                    double score = lineData3TBHN[5] as double
                    double pvalue = lineData3TBHN[7] as double
                    String motifSequence = lineData3TBHN[8]
                    int motifLength = lineData3TBHN[9] as int
                    int locationStart = (lineData3TBHN[11] as int) < 1 ? 1 : (lineData3TBHN[11] as int)
                    // A starting position that is more than 5 figures merges into the high column
                    String locationStartStr = locationStart.toString()
                    if (locationStartStr.length() >= 6) {
                        locationStartStr = locationStartStr.substring(0, 5)
                    }
                    locationStart = locationStartStr as int
                    int locationEnd = locationStart + motifLength - 1
                    if (motifSequence.endsWith("#")) { // it overhangs the protein seq so adjust locationEnd
                        int motifSeqLength = motifSequence.length()
                        int indexCheck = motifSeqLength - 1
                        while (motifSequence[indexCheck] == "#") {
                            indexCheck -= 1
                        }
                        locationEnd = locationEnd - (motifSeqLength - indexCheck) + 1
                    }

                    Prints matchLocation = Prints.buildLocationMatch(
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
    Map<String, Map<String, Match>> matches = [:]
    rawMatches.each { proteinAccession, proteinMatches ->
        List<Prints> sortedMatches = sortMatches(proteinMatches)
        List<Prints> filteredMatches = []
        String currentModelAcc = null
        List<Prints> motifMatchesForCurrentModel = []
        boolean currentMatchesPass = true
        boolean passed = false
        HierarchyEntry currentHierarchyEntry = null
        Set<String> hierarchyEntryIdLimitation = hierarchyMap.keySet() // initialise with all models

        for (Prints rawMatch in sortedMatches) {
            if (currentModelAcc == null || currentModelAcc != rawMatch.modelName) {
                // just started or moved onto a match for a different model
                if (currentModelAcc != null && currentMatchesPass) {
                    passed = selectMatches(
                        motifMatchesForCurrentModel, currentModelAcc,
                        currentHierarchyEntry, hierarchyEntryIdLimitation
                    )
                    if (passed) {
                        filteredMatches.addAll(motifMatchesForCurrentModel)
                        if (currentHierarchyEntry.siblingsIds.length < hierarchyEntryIdLimitation.size()) {
                            hierarchyEntryIdLimitation = new HashSet<>(Arrays.asList(currentHierarchyEntry.siblingsIds))
                        }
                    }
                }

                // reset the values
                currentMatchesPass = true
                motifMatchesForCurrentModel = []  as List<Prints>
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
            for (Prints filteredMatch: filteredMatches) {
                // the modelName has been used up to this point, but we need to convert to the model ID
                Match match = finalMatches.computeIfAbsent(
                    filteredMatch.modelId,
                    {
                        new Match(filteredMatch.modelId, filteredMatch.evalue, filteredMatch.graphscan, new Signature(filteredMatch.modelId, library))
                    }
                )
                Location location = new Location(
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

    def outputFilePath = task.workDir.resolve("prints.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

List<Prints> sortMatches(List<Prints> matches) {
    // This comparator is CRITICAL to the working of PRINTS post-processing
    return matches.sort { matchA, matchB ->
        int evalueComparison = matchA.evalue <=> matchB.evalue
        if (evalueComparison != 0) return evalueComparison

        int modelAccessionComparison = matchA.modelName <=> matchB.modelName
        if (modelAccessionComparison != 0) return modelAccessionComparison

        int motifNumberComparison = matchA.motifNumber <=> matchB.motifNumber
        if (motifNumberComparison != 0) return motifNumberComparison

        int startLocationComparison = matchA.locationStart <=> matchB.locationStart
        if (startLocationComparison != 0) return startLocationComparison

        return matchA.locationEnd <=> matchB.locationEnd
    }
}

boolean selectMatches(  // check if the matches should be selected. Returns a boolean
        List<Match> motifMatchesForCurrentModel,
        String modelName,
        HierarchyEntry hierarchy,
        Set<String> hierarchyEntryIdLimitation
) {
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