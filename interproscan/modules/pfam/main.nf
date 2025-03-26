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
    val datPath

    output:
    tuple val(meta), path("pfam.json")

    exec:
    int MINLENGTH = 8  // minimum length of a fragment
    def outputFilePath = task.workDir.resolve("pfam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "Pfam")
    Map<String, Map<String, Object>> dat = stockholmDatParser(datPath)  // [modelAcc: [clan: str, nested: [str]]]

    def filteredMatches = filterMatches(hmmerMatches, dat)
    def processedMatches = buildFragments(dat, filteredMatches, MINLENGTH)

    json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

def stockholmDatParser(String pfamADatFile) {
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
            parsedDat[accession]?.nested?.add(nestedAcc) ?: (parsedDat[accession] = [nested: [nestedAcc]])
         } else if (line.startsWith("#=GF CL")) {
            String claAcc = line.split()[2]
            parsedDat.computeIfAbsent(accession, { [clan: claAcc] }).clan = claAcc
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

def flatMatchLocations(matches) {
    // Separate locations such that each location is treated as an independent match
    matches.collectMany { modelAccession, match ->
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
}

def filterMatches(Map<String, Map<String, Match>> hmmerMatches, Map<String, Map<String, Object>> dat) {
    /*
        A match is ignored only if it meets all of the following conditions when compared 
        to previously evaluated Pfam matches:
        - It overlaps with another match.
        - Both matches belong to the same clan.
        - Neither match is fully nested within the other.
        Reasoning: Pfam matches aren't supposed to overlap, but families within the 
        same clan are evolutionary related so overlapping is common.
        Therefore, the `nested` field in `pfam_a.dat` is used to whitelist overlaps that 
        can occur between two families.
    */
    Map<String, List<Match>> filteredMatches = [:]
    hmmerMatches.each { seqId, matches ->
         List<Match> allMatches = flatMatchLocations(matches)

        // Sort matches by evalue ASC, score DESC to keep the best matches
        allMatches.sort { a, b ->
            (a.locations[0].evalue <=> b.locations[0].evalue) ?: -(a.locations[0].score <=> b.locations[0].score)
        }

        filteredMatches[seqId] = []
        allMatches.each { match ->
            boolean keep = true
            Map<String, List<String>> candidateMatch = dat[match.modelAccession] ?: [:]
            String candidateMatchClan = candidateMatch?.clan
            filteredMatches[seqId].each { filteredMatch -> // iterate through the matches that have already been chosen
                Map<String, List<String>> filteredMatchInfo = dat[filteredMatch.modelAccession] ?: [:]
                String filteredMatchClan = filteredMatchInfo?.clan
                if (candidateMatchClan == filteredMatchClan) {  // check if both are in the same clan
                    boolean overlapped = isOverlapping(
                        match.locations[0].start, match.locations[0].end,
                        filteredMatch.locations[0].start, filteredMatch.locations[0].end
                    )
                    if (overlapped) {
                        List<String> candidateNested = candidateMatch?.nested ?: []
                        List<String> filteredNested = filteredMatchInfo?.nested ?: []
                        // check if candidates are NOT nested
                        if (!(candidateNested.contains(filteredMatch.modelAccession) || filteredNested.contains(match.modelAccession))) {
                            keep = false
                        }
                    }
                }
            }
            if (keep) {
                filteredMatches[seqId] << match
            }
        }
    }
    return filteredMatches
}

def buildFragments(Map<String, Map<String, Object>> dat,
                               Map<String, Map<String, Match>> filteredMatches,
                               int MINLENGTH) {
    /* Add fragmentLocation objects to the Matches and gather Matches for the same
    model accession into a single Match object */
    Map<String, Map<String, Match>> processedMatches = [:]
    filteredMatches.each { String seqId, List<Match> matches ->
        def matchesAggregated = [:]
        matches.each { Match match ->
            List<String> nestedModels = dat[match.modelAccession]?.nested ?: []
            if (nestedModels) {
                /* Find all matches whose models are listed as nested within the current model 
                in pfam_a.dat AND which overlap the current match.
                Then these matches are converted into a list of maps, that each contain a 
                `start` and `end` field. */
                List<Map<String, Integer>> locationFragments = matches.findAll { otherMatch ->
                    otherMatch.modelAccession in nestedModels &&
                    isOverlapping(otherMatch.locations[0].start, otherMatch.locations[0].end,
                                    match.locations[0].start, match.locations[0].end)
                }.collect { otherMatch ->
                    [start: otherMatch.locations[0].start, end: otherMatch.locations[0].end]
                }
                locationFragments.sort { a, b -> a.start <=> b.start ?: a.end <=> b.end }

                /* Process fragments and identify discontinuous matches.
                The fragmentDcStatus var tracks whether a match is continuous or has discontinuities
                due to nested fragments */
                String fragmentDcStatus = "CONTINUOUS"
                def discontinuousMatchesList = [match]
                locationFragments.each { fragment ->
                    // As we iterate over each fragment, we will attempt to split or adjust matches
                    def newMatchesFromFragment = []
                    discontinuousMatchesList.each { rawDiscontinuousMatch ->
                        List<LocationFragment> fragments = []
                        int newLocationStart = rawDiscontinuousMatch.locations[0].start
                        int newLocationEnd = rawDiscontinuousMatch.locations[0].end
                        int finalLocationEnd = rawDiscontinuousMatch.locations[0].end

                        // The fragment does NOT overlap the match, so the match is left unchanged
                        if (! isOverlapping(newLocationStart, newLocationEnd, fragment['start'], fragment['end'])) {
                            newMatchesFromFragment << rawDiscontinuousMatch
                            return
                        }

                        if (fragment['start'] <= newLocationStart && fragment['end'] >= newLocationEnd) {
                            fragmentDcStatus = "NC_TERMINAL_DISC"
                            fragments.add(new LocationFragment(newLocationStart, newLocationEnd, fragmentDcStatus))
                            rawDiscontinuousMatch.locations[0].fragments = fragments
                            newMatchesFromFragment << rawDiscontinuousMatch
                            return
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

                        if (newLocationEnd - newLocationStart + 1 >= MINLENGTH) {
                            fragments.add(new LocationFragment(newLocationStart, newLocationEnd, fragmentDcStatus))
                            rawDiscontinuousMatch.locations[0].fragments = fragments
                        }
                        newLocationStart = fragment.end + 1
                        if (twoActualRegions) {
                            //deal with final region
                            fragmentDcStatus = "N_TERMINAL_DISC"
                            rawDiscontinuousMatch.locations[0].end = finalLocationEnd  // ensure the 2nd frag extends to the original match's end position
                            // Keep the 2nd region only if it meets the MINLENGTH and store as a new LocationFragment inside discontinuousMatch
                            if (finalLocationEnd - newLocationStart + 1 >= MINLENGTH) {
                                fragments.add(new LocationFragment(newLocationStart, finalLocationEnd, fragmentDcStatus))
                                rawDiscontinuousMatch.locations[0].fragments = fragments
                            }
                        }
                        
                        // If only one valid fragment remains after filtering, update the match start and end to reflect the new fragment boundaries
                        if (rawDiscontinuousMatch.locations[0].fragments.size() == 1) {
                            rawDiscontinuousMatch.locations[0].start = rawDiscontinuousMatch.locations[0].fragments[0].start
                            rawDiscontinuousMatch.locations[0].end = rawDiscontinuousMatch.locations[0].fragments[0].end
                        }
                        
                        newMatchesFromFragment << rawDiscontinuousMatch
                    }
                    // filter out fragment matches that are shorter than MINLENGTH
                    discontinuousMatchesList = newMatchesFromFragment.findAll { it.locations[0].end - it.locations[0].start + 1 >= MINLENGTH }
                }
                discontinuousMatchesList.each { rawMatch ->
                    if (matchesAggregated.containsKey(rawMatch.modelAccession)) {
                        matchesAggregated[rawMatch.modelAccession].locations << rawMatch.locations[0]
                    } else {
                        matchesAggregated[rawMatch.modelAccession] = rawMatch
                    }
                }
            } else if (match.locations[0].end - match.locations[0].start + 1 >= MINLENGTH) {  // filter out matches that are shorter than MINLENGTH
                matchesAggregated.computeIfAbsent(match.modelAccession, { match }).locations << match.locations[0]
            }
        }
        processedMatches[seqId] = matchesAggregated
    }
    return processedMatches
}