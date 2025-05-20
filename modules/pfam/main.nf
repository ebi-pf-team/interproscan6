import groovy.json.JsonOutput

process PARSE_PFAM {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmsearch_out)
    val dirpath
    val datfile

    output:
    tuple val(meta), path("pfam.json")

    exec:
    int MINLENGTH = 8  // minimum length of a fragment
    def outputFilePath = task.workDir.resolve("pfam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "Pfam")
    String datPath = "${dirpath.toString()}/${datfile}"
    Map<String, Map<String, Object>> dat = stockholmDatParser(datPath)  // [modelAcc: [clan: str, nested: [str]]]

    Map<String, List<Match>> filteredMatches = filterMatches(hmmerMatches, dat)
    Map<String, Map<String, Match>> processedMatches = buildFragments(filteredMatches, dat, MINLENGTH)

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

def filterMatches(Map<String, Map<String, Match>> hmmerMatches, Map<String, Map<String, Object>> dat) {
    /*
        Remove overlapping matches, only keeping the one with the best e-value or score.
        Pfam matches are not supposed to overlap, but some families can be evolutionary related,
        especially those belonging to the same clan, so curators can whitelist overlaps using the
        NE (nested) flag.
    */
    Map<String, List<Match>> filteredMatches = [:]
    hmmerMatches.each { seqId, matches ->
        List<Match> allMatches = flattenMatchLocations(matches)

        // Sort matches by evalue ASC, score DESC to keep the best matches
        allMatches.sort { a, b ->
            (a.locations[0].evalue <=> b.locations[0].evalue) ?: -(a.locations[0].score <=> b.locations[0].score)
        }

        filteredMatches[seqId] = []
        allMatches.each { match ->
            boolean keep = true

            Map<String, List<String>> candidateFamily = dat[match.modelAccession] ?: [:]
            List<String> candidateNested = candidateFamily?.nested ?: []

            // Compare the current match with matches already accepted
            filteredMatches[seqId].each { filteredMatch ->
                boolean overlapped = isOverlapping(
                        match.locations[0].start, 
                        match.locations[0].end,
                        filteredMatch.locations[0].start, 
                        filteredMatch.locations[0].end
                )

                if (overlapped) {
                    // Matches are overlapping

                    Map<String, List<String>> filteredFamily = dat[filteredMatch.modelAccession] ?: [:]
                    List<String> filteredNested = filteredFamily?.nested ?: []

                    boolean canBeNested = (candidateNested.contains(filteredMatch.modelAccession)
                                           || filteredNested.contains(match.modelAccession))
                    boolean fullyEnclosed = isFullyEnclosed(
                            match.locations[0].start,
                            match.locations[0].end,
                            filteredMatch.locations[0].start,
                            filteredMatch.locations[0].end
                    )

                    if (!canBeNested || !fullyEnclosed) {
                        keep = false
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

def flattenMatchLocations(matches) {
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

def isOverlapping(location1Start, location1End, location2Start, location2End) {
    Math.max(location1Start, location2Start) <= Math.min(location1End, location2End)
}

def isFullyEnclosed(location1Start, location1End, location2Start, location2End) {
    return (location1Start <= location2Start && location1End >= location2End) 
           || (location2Start <= location1Start && location2End >= location1End)
}

def buildFragments(Map<String, Map<String, Match>> filteredMatches,
                   Map<String, Map<String, Object>> dat,
                   int MINLENGTH) {
    Map<String, Map<String, Match>> processedMatches = [:]
    filteredMatches.each { String seqId, List<Match> matches ->
        def aggregatedMatches = [:]
        matches.each { Match match ->
            // We flattened matches and locations so one location only per match here
            def location = match.locations[0]

            List<String> nestedModels = dat[match.modelAccession]?.nested ?: []
            if (nestedModels) {
                /*
                    Find all locations belonging to matches that:
                    - overlap the current location
                    - belong to models that can be nested within the model of the current match
                */
                List<Map<String, Integer>> overlappingLocations = matches
                    .findAll { otherMatch ->
                        otherMatch.modelAccession in nestedModels 
                        && isFullyEnclosed(
                            otherMatch.locations[0].start, 
                            otherMatch.locations[0].end,
                            location.start, 
                            location.end
                        )
                    }
                    .collect { otherMatch ->
                        [start: otherMatch.locations[0].start, end: otherMatch.locations[0].end]
                    }

                // Sort these locations by position
                overlappingLocations.sort { a, b -> a.start <=> b.start ?: a.end <=> b.end }

                def fragments = []
                int start = location.start

                overlappingLocations.each { otherLocation ->
                    if (otherLocation.start > start && otherLocation.end < location.end) {
                        String status = fragments.isEmpty() ? "C_TERMINAL_DISC" : "NC_TERMINAL_DISC"
                        fragments.add(new LocationFragment(start, otherLocation.start - 1, status))
                        start = otherLocation.end + 1
                    }
                }

                fragments.add(new LocationFragment(start, location.end, "N_TERMINAL_DISC"))
                location.fragments = fragments
            }
            
            if (location.end - location.start + 1 >= MINLENGTH) {
                if (aggregatedMatches.containsKey(match.modelAccession)) {
                    aggregatedMatches[match.modelAccession].locations << location
                } else {
                    aggregatedMatches[match.modelAccession] = match
                }
            }
        }
        processedMatches[seqId] = aggregatedMatches
    }
    return processedMatches
}

def storeMatches(Map<String, Match> aggregatedMatches, List<Match> matches) {
    // Aggregate processed matches under their respective model accessions
    matches.each { match ->
        if (aggregatedMatches.containsKey(match.modelAccession)) {
            aggregatedMatches[match.modelAccession].locations << match.locations[0]
        } else {
            aggregatedMatches[match.modelAccession] = match
        }
    }
    return aggregatedMatches
}