import groovy.json.JsonOutput

process SEARCH_PFAM {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 61295632 --cut_ga \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process PARSE_PFAM {
    label 'run_locally'

    input:
    tuple val(meta), val(hmmsearch_out)
    val seedPath
    val clanPath
    val datPath

    output:
    tuple val(meta), path("pfam.json")

    exec:
    def outputFilePath = task.workDir.resolve("pfam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "Pfam")
    Map<String, Map<String, Object>> dat = stockholmDatParser(datPath)  // [modelAcc: [clan: str, nested: [str]]]
    Map<String, Map<String, Match>> filteredMatches = [:]
    minLength = 8

    // filter matches
    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        // sorting matches by evalue ASC, score DESC to keep the best matches
         // Separate locations such that each location is treated as an independent match
         def allMatches = matches.collectMany { modelAccession, match ->
            modelAccession = modelAccession.split("\\.")[0]
            match.locations.collect { location ->
                Match matchInfo = new Match(
                                modelAccession,
                                match.evalue,
                                match.score,
                                match.bias,
                                match.signature
                            )
                matchInfo.locations = [location]
                matchInfo.signature.accession = modelAccession
                return matchInfo
            }
        }
        // Sort matches by evalue ASC, score DESC
        allMatches.sort { a, b ->
            (a.locations[0].evalue <=> b.locations[0].evalue) ?: -(a.locations[0].score <=> b.locations[0].score)
        }

        filteredMatches[seqId] = []
        allMatches.each { match ->
            boolean keep = true
            Map<String, List<String>> candidateMatch = dat[match.modelAccession] ?: [:]
            String candidateClan = candidateMatch?.["clan"]
            filteredMatches[seqId].each { filteredMatch -> // iterates through the matches already chosen
                Map<String, List<String>> filteredMatchInfo = dat[filteredMatch.modelAccession] ?: [:]
                String filteredMatchClan = filteredMatchInfo?.["clan"]
                if (candidateClan == filteredMatchClan) {  // check if both are on the same clan
                    boolean matchOverlapFiltered = isOverlapping(
                        match.locations[0].start, match.locations[0].end,
                        filteredMatch.locations[0].start, filteredMatch.locations[0].end
                    )
                    if (matchOverlapFiltered) {
                        List<String> candidateNested = candidateMatch?.["nested"] ?: []
                        List<String> filteredNested = filteredMatchInfo?.["nested"] ?: []
                        boolean matchesAreNested = (candidateNested.contains(filteredMatch.modelAccession) || filteredNested.contains(match.modelAccession))
                        if (!matchesAreNested) {
                            keep = false
                        }
                    }
                }
            }
            if (keep) {
                filteredMatches[seqId] << match
            }
        }
        [seqId: filteredMatches[seqId]]
    }

    // build fragments
    processedMatches = filteredMatches.collectEntries { seqId, matches ->
        def matchesAggregated = [:]
        matches.each { match ->
            List<String> nestedModels = dat.get(match.modelAccession)?.nested ?: []
            if (nestedModels) {
                List<Map<String, Integer>> locationFragments = matches.findAll { otherMatch ->
                    otherMatch.modelAccession in nestedModels &&
                    isOverlapping(otherMatch.locations[0].start, otherMatch.locations[0].end,
                                    match.locations[0].start, match.locations[0].end)
                }.collect { otherMatch ->
                    [start: otherMatch.locations[0].start, end: otherMatch.locations[0].end]
                }
                locationFragments.sort { a, b -> a.start <=> b.start ?: a.end <=> b.end }

                String fragmentDcStatus = "CONTINUOUS"
                def rawDiscontinuousMatches = []
                rawDiscontinuousMatches << match
                locationFragments.each { fragment ->
                    def newMatchesFromFragment = []
                    rawDiscontinuousMatches.each { rawDiscontinuousMatch ->
                        List<LocationFragment> fragments = []
                        int newLocationStart = rawDiscontinuousMatch.locations[0].start
                        int newLocationEnd = rawDiscontinuousMatch.locations[0].end
                        int finalLocationEnd = rawDiscontinuousMatch.locations[0].end

                        if (! isOverlapping(newLocationStart, newLocationEnd, fragment['start'], fragment['end'])) {
                            newMatchesFromFragment << rawDiscontinuousMatch
                            return
                        }

                        if (fragment['start'] <= newLocationStart && fragment['end'] >= newLocationEnd) {
                            fragmentDcStatus = "NC_TERMINAL_DISC"
                            fragments.add(new LocationFragment(rawDiscontinuousMatch.locations[0].start, rawDiscontinuousMatch.locations[0].end, fragmentDcStatus))
                            rawDiscontinuousMatch.locations[0].fragments = fragments
                            newMatchesFromFragment << rawDiscontinuousMatch
                            return
                        }

                        if (fragmentDcStatus == "CONTINUOUS") {
                            fragmentDcStatus = null
                        }

                        boolean areSeparateFrags = false
                        if (fragment['start'] <= newLocationStart) {
                            newLocationStart = fragment['end'] + 1
                            fragmentDcStatus = "N_TERMINAL_DISC"
                        } else if (fragment['end'] >= newLocationEnd) {
                            newLocationEnd = fragment['start'] - 1
                            fragmentDcStatus = "C_TERMINAL_DISC"
                        } else if (fragment['start'] > newLocationStart && fragment['end'] < newLocationEnd) {
                            newLocationEnd = fragment['start'] - 1
                            twoActualRegions = true
                            fragmentDcStatus = "C_TERMINAL_DISC"
                        }

                        if (newLocationEnd - newLocationStart + 1 >= minLength) {
                            fragments.add(new LocationFragment(newLocationStart, newLocationEnd, fragmentDcStatus))
                            rawDiscontinuousMatch.locations[0].fragments = fragments
                        }
                        newLocationStart = fragment.end + 1
                        if (twoActualRegions) {
                            //deal with final region
                            fragmentDcStatus = "N_TERMINAL_DISC"
                            rawDiscontinuousMatch.locations[0].end = finalLocationEnd
                            if (finalLocationEnd - newLocationStart + 1 >= minLength) {
                                fragments.add(new LocationFragment(newLocationStart, finalLocationEnd, fragmentDcStatus))
                                rawDiscontinuousMatch.locations[0].fragments = fragments
                            }
                        }
                        if (rawDiscontinuousMatch.locations[0].fragments.size() == 1) {
                            rawDiscontinuousMatch.locations[0].start = rawDiscontinuousMatch.locations[0].fragments[0].start
                            rawDiscontinuousMatch.locations[0].end = rawDiscontinuousMatch.locations[0].fragments[0].end
                        }
                        newMatchesFromFragment << rawDiscontinuousMatch
                    }
                    // filter out fragment matches that are shorter than minLength
                    rawDiscontinuousMatches = newMatchesFromFragment.findAll { it.locations[0].end - it.locations[0].start + 1 >= minLength }
                }
                rawDiscontinuousMatches.each { rawMatch ->
                    if (matchesAggregated.containsKey(rawMatch.modelAccession)) {
                        matchesAggregated[rawMatch.modelAccession].locations << rawMatch.locations[0]
                    } else {
                        matchesAggregated[rawMatch.modelAccession] = rawMatch
                    }
                }
            } else if (match.locations[0].end - match.locations[0].start + 1 >= minLength) {  // filter out matches that are shorter than minLength
                if (matchesAggregated.containsKey(match.modelAccession)) {
                    matchesAggregated[match.modelAccession].locations << match.locations[0]
                } else {
                    matchesAggregated[match.modelAccession] = match
                }
            }
        }

        [seqId: matchesAggregated]
    }

    json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

Map<String, Map<String, Object>> stockholmDatParser(String pfamADatFile) {
    /* Retrieve nested models and clan classifications.
    E.g. [ PF00026:[nested:[PF03489, PF05184], clan:CL0129], PF06826:[clan:CL0064, nested:[]] ]
    */
    Map<String, String> name2acc = [:]
    Map<String, List<String>> parsedDat = [:]
    String name = null
    String accession = null
    new File(pfamADatFile).eachLine { rawLine ->
        def line = decode(rawLine.bytes)
        if (line.startsWith("#=GF ID")) {
            name = line.split()[2]
        } else if (line.startsWith("#=GF AC")) {
            accession = line.split()[2].split("\\.")[0]
            name2acc[name] = accession
        } else if (line.startsWith("#=GF NE")) {
            String nestedAcc = line.split()[2]
            parsedDat[accession]?.nested?.add(nestedAcc) ?:
                (parsedDat[accession] = [nested: ([nestedAcc])])
         } else if (line.startsWith("#=GF CL")) {
            String claAcc = line.split()[2]
            if (parsedDat[accession]) {
                parsedDat[accession].clan = claAcc
            } else {
                parsedDat[accession] = [clan: claAcc]
            }
        }
    }
    // convert nested 'names' to 'acc'
    parsedDat.each { acc, info ->
        def nestedNames = info.nested ?: []
        def nestedAccessions = nestedNames.collect { name2acc[it] }.findAll { it != null }
        info.nested = nestedAccessions.unique()
    }

    return parsedDat
}

def decode(byte[] b) {
    try {
        return new String(b, "UTF-8").trim()
    } catch (Exception e) {
        return new String(b, "ISO-8859-1").trim()
    }
}

def isOverlapping(location1Start, location1End, location2Start, location2End) {
    Math.max(location1Start, location2Start) <= Math.min(location1End, location2End)
}