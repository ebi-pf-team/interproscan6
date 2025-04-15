import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process REPRESENTATIVE_LOCATIONS {
    executor 'local'

    input:
    tuple val(meta), val(matchesPath)

    output:
    tuple val(meta), path("matches_with_representative.json")

    exec:
    int MAX_DOMS_PER_GROUP = 20 // only consider N "best" locations otherwise there are too many comparisons (2^locs)
    float DOM_OVERLAP_THRESHOLD = 0.3
    List<String> REPR_TYPE = ["family", "domain"]

    def matchesMap = new ObjectMapper().readValue(new File(matchesPath.toString()), Map)
    matchesMap.each { String md5, Map matchesInMap ->
        // Serialise the matches so we don't need to edit the map manually later
        Map<String, Match> currentMatches = [:]
        matchesInMap.each { String modelAcc, Map matchMap ->
            currentMatches[modelAcc] = Match.fromMap(matchMap)
        }

        // Look for representatives for matches of a specific type, e.g. "Domain":
        REPR_TYPE.each { String reprType ->
            // Gather relevant locations. Only matches from the relevant dbs will have a representativeInfo object
            List<CandidateLocation> candidateLocations = getCandidateLocations(currentMatches, reprType)

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
                    }.take(MAX_DOMS_PER_GROUP)

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
                                if (locationsOverlap(loc1.residues, grp[j].residues, DOM_OVERLAP_THRESHOLD)) {
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
    }

    def outputFilePath = task.workDir.resolve("matches_with_representative.json")
    def json = JsonOutput.toJson(matchesMap)
    new File(outputFilePath.toString()).write(json)
}

List<CandidateLocation> getCandidateLocations(Map matches, String reprType) {
    List<CandidateLocation> candidateLocations = []
    matches.each { String modelAccession, Match match ->
        if (match.representativeInfo?.type == reprType) {
            match.locations.each { Location loc ->
                CandidateLocation candidate = new CandidateLocation(loc, match.representativeInfo.rank)
                candidateLocations.add(candidate)
            }
        }
    }
    return candidateLocations
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
