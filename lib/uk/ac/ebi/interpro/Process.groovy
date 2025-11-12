package uk.ac.ebi.interpro

import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper
import uk.ac.ebi.interpro.CandidateLocation
import uk.ac.ebi.interpro.Entry
import uk.ac.ebi.interpro.GoXRefs
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.PathwayXRefs
import uk.ac.ebi.interpro.RepresentativeInfo


class Process {
    static void combineMatches(List<String> inputPaths, String outputPath) {
        def mapper = new ObjectMapper()
        def aggregatedMatches = [:]  // seqMd5 -> (modelAcc -> match)

        inputPaths.each { path ->
            def file = new File(path)
            def data = mapper.readValue(file, Map)
            data.each { seqMd5, matches ->
                aggregatedMatches.computeIfAbsent(seqMd5) { [:] }.putAll(matches)
            }
        }

        def jf = new JsonFactory()
        new File(outputPath).withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, aggregatedMatches)
            gen.close()
        }
    }

    static void addXrefs(String inputPath, Map dbReleases, boolean addGoTerms, 
                         boolean addPathways, String pantherPaintDirectory, String outputPath) {
        ObjectMapper mapper = new ObjectMapper()
        def (databaseInfo, entries, ipr2go, goInfo, ipr2pa, paInfo) = [null, null, null, null, null, null]
        if (dbReleases?.interpro?.dirpath != null) {
            String interproDir = dbReleases.interpro.dirpath.toString()
            String databasesPath = "${interproDir}/databases.json"
            File databasesJson = new File(databasesPath)
            databaseInfo = mapper.readValue(databasesJson, Map)
            String entriesPath = "${interproDir}/entries.json"
            File entriesJson = new File(entriesPath)
            entries = mapper.readValue(entriesJson, Map)
            (ipr2go, goInfo) = loadXRefFiles("${interproDir}/goterms")
            if (addPathways) {
                (ipr2pa, paInfo) = loadXRefFiles("${interproDir}/pathways")
            }
        }

        def matches = mapper.readValue(new File(inputPath), Map.class)
        matches.each { String seqMd5, Map seqMatches ->
            seqMatches.each { modelAcc, matchMap ->
                Match match = Match.fromMap(matchMap)  // convert Map to Match object

                if (match.source != "InterPro-N" && entries) {
                    String signatureAcc = match.signature.accession
                    def signatureInfo = entries[signatureAcc] ?: entries[modelAcc]

                    // Update library version
                    def version = (match.signature.signatureLibraryRelease.version == "null") ? null : match.signature.signatureLibraryRelease.version
                    if (!version && signatureInfo != null) {
                        match.signature.signatureLibraryRelease.version = databaseInfo[signatureInfo["database"]]
                    } else if (match.signature.signatureLibraryRelease.library == "PIRSR") {
                        match.signature.signatureLibraryRelease.version = databaseInfo["PIRSR"]
                    }

                    // Handle PANTHER data
                    if (match.signature.signatureLibraryRelease.library == "PANTHER") {
                        updatePantherData(match, dbReleases.panther.dirpath, pantherPaintDirectory, signatureAcc, entries, addGoTerms ? goInfo : null)
                    }

                    // Update signature info
                    if (signatureInfo != null) {
                        match.signature.name = signatureInfo["name"]
                        match.signature.description = signatureInfo["description"]
                        String sigType = signatureInfo["type"]
                        match.signature.setType(sigType)

                        if (signatureInfo["representative"] != null) {
                            match.representativeInfo = new RepresentativeInfo(
                                signatureInfo["representative"]["type"],
                                signatureInfo["representative"]["index"]
                            )
                        }

                        // Handle InterPro data
                        String interproAcc = signatureInfo["integrated"]
                        if (interproAcc != null) {
                            def entryInfo = entries[interproAcc]
                            assert entryInfo != null
                            match.signature.entry = new Entry(
                                interproAcc, entryInfo["name"], entryInfo["description"], entryInfo["type"]
                            )
                            addXRefs(match, interproAcc, addGoTerms ? ipr2go : null, goInfo, ipr2pa, paInfo)
                        }
                    }

                    // update the match
                    seqMatches[modelAcc] = match
                }  // end of if (entries)
            } // end of seqMatches
        }  // end of Json reader / seq Id

        def jf = new JsonFactory()
        new File(outputPath).withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, matches)
            gen.close()
        }
    }

    static def loadXRefFiles(String prefix) {
        ObjectMapper mapper = new ObjectMapper()
        File iprFile = new File("${prefix}.ipr.json")
        File infoFile = new File("${prefix}.json")

        Map iprData = mapper.readValue(iprFile, Map)
        Map infoData = mapper.readValue(infoFile, Map)

        return [iprData, infoData]
    }

    static void updatePantherData(Match match, String pantherDir, String paintAnnoDir, 
                                 String signatureAcc, Map entries, Map goInfo) {
        Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
        File paintAnnotationFile = new File("${pantherDir.toString()}/${paintAnnoDir}/${signatureAcc}.json")
        assert paintAnnotationFile.exists()
        def paintAnnotationsContent = new ObjectMapper().readValue(paintAnnotationFile, Map)
        String nodeId = match.treegrafter.ancestralNodeID
        def nodeData = paintAnnotationsContent[nodeId]
        if (nodeData != null) {
            match.treegrafter.subfamilyAccession = (nodeData[0] == null) ? null : nodeData[0]
            if (nodeData[1] != null && goInfo?.terms != null) {
                nodeData[1].split(',').each { goTermId ->
                    def term = goInfo.terms[goTermId]
                    if (term) {
                        def (goTermName, goTermType) = term
                        match.treegrafter.addGoXRefs(
                            new GoXRefs(
                                goTermName,
                                "GO",
                                GO_PATTERN[goTermType],
                                goTermId
                            )
                        )
                    }
                }
            }
            match.treegrafter.proteinClass = (nodeData[2] == null) ? null : nodeData[2]
            match.treegrafter.graftPoint = (nodeData[3] == null) ? null : nodeData[3]
        }

        if (entries[signatureAcc]) {
            match.treegrafter.subfamilyName = entries[signatureAcc]["name"]
            match.treegrafter.subfamilyDescription = entries[signatureAcc]["description"]
        }
    }

    static void addXRefs(Match match, String interproAcc, def ipr2go, def goInfo, def ipr2pa, def paInfo) {
        Map<String,String> GO_PATTERN = ["P": "BIOLOGICAL_PROCESS", "C": "CELLULAR_COMPONENT", "F": "MOLECULAR_FUNCTION"]
        Map<String,String> PA_PATTERN = ["t": "MetaCyc", "w": "UniPathway", "k": "KEGG", "r": "Reactome"]
        if (ipr2go != null && goInfo != null && ipr2go[interproAcc] != null) {
            ipr2go[interproAcc].each { goId ->
                goId = goId
                match.signature.entry.addGoXRefs(
                    new GoXRefs(goInfo["terms"][goId][0], "GO", GO_PATTERN[goInfo["terms"][goId][1]], goId)
                )
            }
        }
        if (ipr2pa != null && paInfo != null && ipr2pa[interproAcc] != null) {
            ipr2pa[interproAcc].each { paId ->
                paId = paId
                match.signature.entry.addPathwayXRefs(
                    new PathwayXRefs(paInfo[paId][1], PA_PATTERN[paInfo[paId][0]], paId)
                )
            }
        }
    }

    static void selectRepresentativeLocations(String inputPath, String outputPath) {
        // only consider N "best" locations otherwise there are too many comparisons (2^locs)
        int maxLocationsPerGroup = 20
        float locationOverlapThreshold = 0.3
        List<String> representativeTypes = ["family", "domain"]

        ObjectMapper mapper = new ObjectMapper()
        def matches = mapper.readValue(new File(inputPath), Map.class)
        matches.each { String seqMd5, Map rawMatches ->
            // Deserialise all matches
            Map<String, Match> seqMatches = [:]
            rawMatches.each { String modelAcc, Map rawMatch ->
                seqMatches[modelAcc] = Match.fromMap(rawMatch)
            }

            // Look for representatives for matches of a specific type, e.g. "Domain":
            representativeTypes.each { String reprType ->
                // Gather relevant locations. Only matches from the relevant dbs will have a representativeInfo object
                List<CandidateLocation> candidateLocations = getCandidateLocations(seqMatches, reprType)

                // Identify/select representative domains
                if (!candidateLocations.isEmpty()) {
                    // Sort based on location position
                    candidateLocations.sort { CandidateLocation loc1, CandidateLocation loc2 ->
                        int delta = loc1.location.start - loc2.location.start
                        delta != 0 ? delta : loc1.location.end - loc2.location.end
                    }

                    // Group domains together
                    List<Location> groups = new ArrayList<>()
                    List<Location> group = new ArrayList<>()
                    group.add(candidateLocations[0])
                    int stop = candidateLocations[0].location.end
                    if (candidateLocations.size() > 1) {
                        for (CandidateLocation candidate : candidateLocations[1..-1]) {
                            if (candidate.location.start <= stop) {
                                group.add(candidate)
                                stop = Math.max(stop, candidate.location.end)
                            } else {
                                groups.add(group)
                                group = [candidate]
                                stop = candidate.location.end
                            }
                        }
                    }
                    groups.add(group) // Add the last group

                    // Select representative domain in each group
                    groups.each { grp ->
                        grp = grp.sort { a, b ->
                            // compare the number of residues covered by each match
                            int resComparison = b.residues.size() - a.residues.size()
                            // if their coverage is the same, use the database rank
                            resComparison != 0 ? resComparison : a.representativeRank - b.representativeRank
                        }.take(maxLocationsPerGroup)

                        // Process representative domains in the group
                        if (grp.size() > 1) {
                            Map<Integer, Set<Integer>> graph = (0..<grp.size()).collectEntries { i ->
                                /*
                                `(0..<x) - y` creates a range from 0 (inclusive) to x (exclusive)
                                `- i` removes the value represented by `i` in the range
                                so `(0..<5) - 2`, the result is `[0, 1, 3, 4]`
                                */
                                [i, new HashSet<>((0..<grp.size()) - i)]
                            }
                            grp.eachWithIndex { loc1, i ->
                                (i + 1..<grp.size()).each { j ->
                                    if (locationsOverlap(loc1.residues, grp[j].residues, locationOverlapThreshold)) {
                                        graph[i].remove(Integer.valueOf(j))
                                        graph[j].remove(Integer.valueOf(i))
                                    }
                                }
                            }

                            // Select the best subgroup
                            Set<Set<Integer>> subgroups = getValidSets(graph)
                            def bestSubgroup = null
                            int maxCoverage = 0
                            int maxPfams = 0
                            subgroups.each { subgrp ->
                                Set coverage = new HashSet<>()
                                int pfams = 0
                                List currentGrp = []
                                subgrp.each { i ->
                                    def loc = grp[i]
                                    coverage.addAll(loc.residues)
                                    if (loc.representativeRank == 0) { pfams++ }
                                    currentGrp.add(loc)
                                }
                                int currentCoverage = coverage.size()
                                if (currentCoverage > maxCoverage || (currentCoverage == maxCoverage && pfams > maxPfams)) {
                                    maxCoverage = currentCoverage
                                    maxPfams = pfams
                                    bestSubgroup = currentGrp
                                }
                            }
                            bestSubgroup.each { candidate -> candidate.location.representative = true }
                        } else {
                            grp.each { candidate -> candidate.location.representative = true }
                        }
                    }
                }
            }

            matches[seqMd5] = seqMatches
        }

        def jf = new JsonFactory()
        new File(outputPath).withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, matches)
            gen.close()
        }
    }

    static List<CandidateLocation> getCandidateLocations(Map matches, String reprType) {
        List<CandidateLocation> candidateLocations = []
        matches.each { String modelAccession, Match match ->
            if (match.source != "InterPro-N" && match.representativeInfo?.type == reprType) {
                match.locations.each { Location loc ->
                    CandidateLocation candidate = new CandidateLocation(loc, match.representativeInfo.rank)
                    candidateLocations.add(candidate)
                }
            }
        }
        return candidateLocations
    }

    static boolean locationsOverlap(Set<Integer> loc1Residues, Set<Integer> loc2Residues, float threshold) {
        int overlap = loc1Residues.intersect(loc2Residues).size()
        return overlap > 0 && (overlap / Math.min(loc1Residues.size(), loc2Residues.size())) >= threshold
    }

    static Set<Set<Integer>> getValidSets(Map<Integer, Set<Integer>> graph) {
        Set<Set<Integer>> allValidSets = new HashSet<>()
        def setIsValid = { candidate ->
            candidate.every { a -> candidate.every { b -> a == b || graph[a].contains(b) } }
        }
        /* Closures must reference variables that are already defined.
        This is because closures capture vars and their scope during definition.
        Defining buildValidSets first to prevent a MissingPropertyException error. */
        def buildValidSets
        buildValidSets = { currentSet, remainingNodes ->
            if (setIsValid(currentSet)) {
                if (!remainingNodes) {
                    allValidSets << currentSet.toSet()
                } else {
                    def currentNode = remainingNodes[0]
                    def restNodes = remainingNodes.size() > 1 ? remainingNodes[1..-1] : []
                    buildValidSets(currentSet + [currentNode], restNodes)
                    buildValidSets(currentSet, restNodes)
                }
            }
        }

        buildValidSets([], graph.keySet().toList())
        return allValidSets
    }
}
