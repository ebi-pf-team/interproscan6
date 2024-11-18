import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process REPRESENTATIVE_DOMAINS {
    label 'analysis_parser'

    input:
    val matchesPath

    output:
    path("matches_repr_domains.json")

    script:
    int MAX_DOMS_PER_GROUP = 20 // only consider N "best" domains otherwise there are too many comparisons (2^domains)
    float DOM_OVERLAP_THRESHOLD = 0.3
    List<String> REPR_DOM_DBS = ["cdd", "ncbifam", "pfam", "prositeprofiles", "smart"]
    def allMatches = []
    JsonSlurper jsonSlurper = new JsonSlurper()
    def matchesMap = jsonSlurper.parse(matchesPath.toFile())
    matchesMap.each { seqData ->  // keys in seqData: sequence, md5, matches, xref, id
        // Gather relevant locations
        List seqDomains = []
        seqData["matches"].each { modelAccession, matchMap ->
            if (matchMap["signature"]["signatureLibraryRelease"]["library"] in REPR_DOM_DBS) {
                Match match = Match.fromMap(matchMap)
                match.locations.each { loc ->
                    loc.representativeRank = match.representativeInfo.rank
                    loc.sortFragments()
                    loc.getResidues()
                    seqDomains.add(loc)
                }
            }
        }

        // Identify/select representative domains
        if (!seqDomains.isEmpty()) {
        for (Location loc: seqDomains) {
            println("loc ${loc} - ${loc.residues} -- ${loc.fragments}")
        }
            // Sort based on location position
            seqDomains.sort { Location loc1, Location loc2 ->
                int delta = loc1.location.start - loc2.location.start
                delta != 0 ? delta : loc1.location.end - loc2.location.end
            }

            // Group domains together
            List groups = []
            List group = [seqDomains[0]]
            int stop = seqDomains[0].fragments[-1].end
            if (seqDomains.size() > 1) {
                for (Location loc : seqDomains[1..-1]) {
                    if (loc.fragments[0].start <= stop) {
                        group.add(loc)
                        stop = Math.max(stop, loc.fragments[-1].end)
                    } else {
                        groups.add(group)
                        group = [loc]
                        stop = loc.fragments[-1].end
                    }
                }
            }
            groups.add(group) // Add the last group

            // Select representative domain in each group
            groups.each { grp ->
                grp = grp.sort { a, b ->
                    int resComparison = b.residues.size() <=> a.residues.size()  // number fragments in desending order
                    resComparison != 0 ? resComp : a.rank <=> b.rank             // rank in ascending order
                }.take(MAX_DOMS_PER_GROUP)

                // Process representative domains in the group
                if (grp.size() > 1) {
                    Map<Integer, Set<Integer>> graph = (0..<grp.size()).collectEntries { i -> [i, (0..<grp.size()) - i] }

                    grp.eachWithIndex { loc1, i ->
                        (i + 1..<grp.size()).each { j ->
                            if (locationsOverlap(loc1.residues, grp[j].residues, DOM_OVERLAP_THRESHOLD)) {
                                println("i ${i} j ${j} graph ${graph}")
                                graph[i].remove(j)
                                graph[j].remove(i)
                            }
                        }
                    }
                    println("------")

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
                    bestSubgroup.each { loc -> loc.representative = true }
                } else {
                    grp.each { loc -> loc.representative = true }
                }
            }
        }
        // Ensure updates are reflected in seqData
        seqData["matches"].each { modelAccession, matchMap ->
            matchMap["locations"].each { loc ->
                def updatedLoc = seqDomains.find { it.start == loc.start && it.end == loc.end }
                if (updatedLoc) { loc.representative = updatedLoc.representative }
            }
        }
        allMatches.add(seqData)
    }
    def outputFilePath = task.workDir.resolve("matches_repr_domains.json")
    def json = JsonOutput.toJson(allMatches)
    new File(outputFilePath.toString()).write(json)
}

boolean locationsOverlap(Set<Integer> loc1Residues, Set<Integer> loc2Residues, float threshold) {
    int overlap = loc1Residues.intersect(loc2Residues).size()
    return overlap > 0 && (overlap / Math.min(loc1Residues.size(), loc2Residues.size())) >= threshold
}

Set<Set<Integer>> getValidSets(Map<Integer, Set<Integer>> graph) {
    Set<Set<Integer>> allValidSets = new HashSet<>()
    def setIsValid = { candidate ->
        candidate.every { a -> candidate.every { b -> a == b || graph[a].contains(b) } }
    }

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
