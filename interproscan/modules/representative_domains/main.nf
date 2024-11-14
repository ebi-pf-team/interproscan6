import groovy.json.JsonOutput

process REPRESENTATIVE_DOMAINS {
    label 'analysis_parser'

    input:
    val matchesPath

    output:
    path("matches_repr_domains.json")

    exec:
    /* We only consider the "best" N domains of the group,
    otherwise the number of possible combinations/sets is too high,
    for example, if M domains, the max number of combinations is 2^M. */
    MAX_DOMS_PER_GROUP = 20
    DOM_OVERLAP_THRESHOLD = 0.3
    REPR_DOM_DBS = ["cdd", "ncbifam", "pfam", "prositeprofiles", "smart"]

    JsonSlurper jsonSlurper = new JsonSlurper()
    JsonSlurper.parse(matchesPath).each { seqData ->
        // keys in seqData: sequence, md5, matches, xref, id
        // Gather relevant locations
        List<Location> seqDomains = []
        seqData["matches"].each { modelAccession, matchMap ->
            if (matchMap["signature"]["signatureLibraryRelease"]["library"] in REPR_DOM_DBS) {
                Match match = Match.fromMap(matchMap)
                for (Location loc: match.locations) {
                    loc.representativeRank = match.representativeInfo.rank
                    seqDomains.add(loc)
                }
            }
        }

        // Identify/select representative domains
        if (!seqDomains.isEmpty()) {
            // sort by boundaries
            seqDomains.sort { loc1, loc2 ->
                def loc1Start = loc1.fragments[0].start
                def loc1End = loc1.fragments[-1].end
                def loc2Start = loc2.fragments[0].start
                def loc2End = loc2.fragments[-1].end
                loc1Start <=> loc2Start ?: loc1End <=> loc2End
            }

            // Group overlapping domains
            int stop = seqDomains[0].fragments[-1].end
            List<List<Location>> groups = []
            List<Location> group = [seqDomains[0]
            for (Location loc: seqDomains[1:]) {
                start = loc.fragments[0].start
                if (start <= stop) {
                    group.add(loc)
                    stop = Math.max(stop, loc.fragments[-1].end)
                } else {
                    groups.add(group)
                    group = [loc]
                    stop = loc.fragments[-1].end
                }
            }
            groups.add(group)

            // Select representative domain in each group
            for (group: groups) {
                List<Location> reprGroup = group.sort { a, b ->
                    int resComp = getResidues(b).size() <=> getResidues(a).size()
                    resComp != 0 ? resComp : a.rank <=> b.rank
                }.take(MAX_DOMS_PER_GROUP)
                Set<int> nodes = (0..<reprGroup.size()).toSet()
                Map<int, Set<int>> graph = nodes.collectEntries { i -> [1, nodes - i]}
                // e.g. nodes = [0,1,2]
                // e.g. graph = [0: [1,2], 1: [0,2], 2: [0,1]]
                reprGroup.eachWithIndex { loc1, i ->
                    (i + 1..<reprGroup.size()).each { j ->
                        def loc2 = reprGroup[j]
                        if (locationsOverlap(getResidues(loc1), getResidues(loc2), DOM_OVERLAP_THRESHOLD)) {
                            graph[i].remove(j)
                            graph[j].remove(i)
                        }
                    }
                }

                Set<Set<int>> subgroups = getValidSets(graph)

                // Find the best combinations
                int maxCoverage = 0
                int maxPfams = 0
                def bestSubgroup = null
                for (Set<int> grp: subgroups) {
                    Set coverage = new HashSet<>()
                    int pfams = 0
                    List subgrp = []
                    subgrp.each { i ->
                        def loc = grp[i]
                        coverage.addAll(getResidues(loc))
                        if (loc.representativeRank == 0) { pfams += 1 }
                        subgrp.add(loc)
                    }
                    int currentCoverage = coverage.size()
                    if (currentCoverage > maxCoverage || pfams > maxPfams) {
                        maxCoverage = currentCoverage
                        maxPfams = pfams
                        bestSubgroup = subgrp
                    }
                }

                for (loc: bestSubgroup) { loc.representative = true }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("matches_repr_domains.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

Set<Integer> getResidues(Location location) {
    // ???
    Set<Integer> residues = new HashSet<>()
    location.fragments.each { frag -> residues.addAll((frag.start..frag.end).toSet()) }
    return residues
}

boolean locationsOverlap(Location loc1Residues, Location loc2Residues, float threshold) {
    overlap = loc1Residues.size() <=> loc1Residues.size()
    return overlap > 0 && (overlap / Math.min(loc1Residues.size(), loc2Residues.size()) >= threshold
}

Set getValidSets(Map<int, Set<int>> graph) {
    def allValidSets = []
    // Iterate over each pair of nodes in candidate to ensure each node is connected to every other node.
    def setIsValid = { candidate -> candidate.every { a -> candidate.every { b -> a == b || graph[b].contains(a) } } }

    def buildValidSets
    buildValidSets = { currentSet, remainingNodes ->
        if (setIsValid(currentSet)) {
            if (!remainingNodes) {
                allValidSets << currentSet.toSet()
            } else {
                def (currentNode, *restNodes) = remainingNodes
                // Consider both including and excluding the currentNode in/from the currentSet
                buildValidSets(currentSet + [currentNode], restNodes)
                buildValidSets(currentSet, restNodes)
            }
        }
    }
    // start with an empty currentSet, and all nodes from the graph
    buildValidSets([], graph.keySet().toList())
    return allValidSets
}
