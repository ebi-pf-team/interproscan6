import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process REPRESENTATIVE_DOMAINS {
    label 'analysis_parser'

    input:
    val matchesPath

    output:
    path "matches_repr_domains.json"

    exec:
    /* We only consider the "best" N domains of the group,
    otherwise the number of possible combinations/sets is too high,
    for example, if M domains, the max number of combinations is 2^M. */
    int MAX_DOMS_PER_GROUP = 20
    float DOM_OVERLAP_THRESHOLD = 0.3
    List<String> REPR_DOM_DBS = ["cdd", "ncbifam", "pfam", "prositeprofiles", "smart"]

    JsonSlurper jsonSlurper = new JsonSlurper()
    jsonSlurper.parse(matchesPath).each { seqData ->
        // keys in seqData: sequence, md5, matches, xref, id
        // Gather relevant locations
        List seqDomains = []
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
            List groups = []
            List group = [seqDomains[0]]
            for (Location loc: seqDomains[1..-1]) {
                int start = loc.fragments[0].start
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
            for (grp in groups) {
                def reprGroup = grp.sort { a, b ->
                    int resComp = getResidues(b).size() <=> getResidues(a).size()
                    resComp != 0 ? resComp : a.rank <=> b.rank
                }.take(MAX_DOMS_PER_GROUP)
                Set<Integer> nodes = (0..<reprGroup.size()).toSet()
                Map<Integer, Set<Integer>> graph = nodes.collectEntries { i -> [i, nodes - i]}
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

                Set<Set<Integer>> subgroups = getValidSets(graph)

                // Find the best combinations
                int maxCoverage = 0
                int maxPfams = 0
                def bestSubgroup = null
                for (Set<Integer> currentGrp: subgroups) {
                    Set coverage = new HashSet<>()
                    int pfams = 0
                    List subgrp = []
                    currentGrp.each { i ->
                        def loc = reprGroup[i]
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

                bestSubgroup.each { loc -> loc.representative = true }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("matches_repr_domains.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

Set<Integer> getResidues(Location location) {
    Set<Integer> residues = new HashSet<>()
    location.fragments.each { frag -> residues.addAll((frag.start..frag.end).toSet()) }
    return residues
}

boolean locationsOverlap(Set<Integer> loc1Residues, Set<Integer> loc2Residues, float threshold) {
    int overlap = loc1Residues.intersect(loc2Residues).size() // Declared 'overlap' variable
    return overlap > 0 && (overlap / Math.min(loc1Residues.size(), loc2Residues.size())) >= threshold // Closed this line
}

Set<Set<Integer>> getValidSets(Map<Integer, Set<Integer>> graph) {
    def allValidSets = []
    // Iterate over each pair of nodes in candidate to ensure each node is connected to every other node.
    def setIsValid = { candidate ->
        candidate.every { a ->
            candidate.every { b ->
                a == b || graph[b].contains(a)
            }
        }
    }

    def buildValidSets = { currentSet, remainingNodes ->
        if (setIsValid(currentSet)) {
            if (!remainingNodes) {
                allValidSets << currentSet.toSet()
            } else {
                def currentNode = remainingNodes[0]
                def restNodes = remainingNodes[1..-1] ?: []
                // Consider both including and excluding the currentNode in/from the currentSet
                buildValidSets(currentSet + [currentNode], restNodes)
                buildValidSets(currentSet, restNodes)
            }
        }
    }

    buildValidSets([], graph.keySet().toList())
    return allValidSets.toSet()
}
