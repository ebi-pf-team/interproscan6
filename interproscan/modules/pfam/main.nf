import groovy.json.JsonOutput

process SEARCH_PFAM {
    label 'hmmer_runner'

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
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmsearch_out)
    val seedPath
    val clanPath
    val datPath

    output:
    tuple val(meta), path("pfam.json")

    exec:
    def outputFilePath = task.workDir.resolve("pfam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString())

    seeds = stockholmSeedParser(seedPath)
    nestedInfo = stockholmClansParser(clanPath, seeds)
    dat = stockholmDatParser(datPath)

    Map<String, Map<String, Match>> filteredMatches = [:]

    // filter matches
    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        // sorting matches by evalue ASC, score DESC to keep the best matches
        Map<String, Match> sortedMatches = matches.sort { a, b ->
            (a.value.evalue <=> b.value.evalue) ?: -(a.value.score <=> b.value.score)
        }
        filteredMatches[seqId] = [:]
        sortedMatches.each { modelAccession, match ->
            modelAccession = modelAccession.split("\\.")[0]
            match.modelAccession = modelAccession
            boolean keep = true
            Map<String, List<String>> candidateMatch = nestedInfo[modelAccession] ?: [:]
            String candidateClan = candidateMatch?.get("clan", null)

            filteredMatches[seqId].each { filteredAcc, filteredMatch -> // iterates through the matches already chosen
                Map<String, List<String>> filteredMatchInfo = nestedInfo[filteredAcc] ?: [:]
                String filteredMatchClan = filteredMatchInfo?.get("clan", null)
                if (candidateClan && candidateClan == filteredMatchClan) {  // check if both are on the same clan
                    List<Location> notOverlappingLocations = match.locations.findAll { candidateLocation ->
                        boolean isOverlapping = filteredMatch.locations.any { filteredLocation ->
                            isOverlapping(candidateLocation, filteredLocation)
                        }
                        if (isOverlapping) {
                            List<String> candidateNested = candidateMatch.get("nested", [])
                            List<String> filteredNested = filteredMatchInfo.get("nested", [])
                            boolean matchesAreNested = candidateNested && filteredNested && candidateNested.intersect(filteredNested)
                            matchesAreNested  // if they are nested, keep both
                        } else {
                            true  // keep if not overlapping
                        }
                    }
                    if (!notOverlappingLocations) {
                        keep = false
                    } else {
                        match.locations = notOverlappingLocations
                    }
                }
            }
            if (keep) {
                filteredMatches[seqId][modelAccession] = match
            }
        }
    }

    // build fragments
    minResiduesLength = 8
    processedMatches = filteredMatches.collectEntries { seqId, matches ->
        [(seqId): matches.findAll { modelAccession, match -> match.locations }
        .collectEntries { modelAccession, match ->
            List<String> nestedModels = dat.get(modelAccession, [])
            if (nestedModels) {
                List<Map<String, Integer>> locationFragments = matches.findAll { otherAccession, _ -> otherAccession in nestedModels }
                    .collectMany { _, otherMatch ->
                        otherMatch.locations.findAll { loc -> isOverlapping(loc, match.locations[0]) }
                                            .collect { loc -> [start: loc.start, end: loc.end] }
                    }
                    .sort { a, b -> a['start'] <=> b['start'] ?: a['end'] <=> b['end'] }
                String fragmentDcStatus = "CONTINUOUS"
                locationFragments.each { fragment ->
                    def (newLocationStart, newLocationEnd, finalLocationEnd, updatedStatus, twoActualRegions) =
                        createFragment(match, fragment, fragmentDcStatus, false)

                    match.locations.each { loc ->
                        loc.fragments = [new LocationFragment(newLocationStart, newLocationEnd, updatedStatus)]
                        if (twoActualRegions) {
                            loc.fragments << new LocationFragment(fragment['end'] + 1, finalLocationEnd, "N_TERMINAL_DISC")
                        }
                    }
                }
            }
            List<Location> validLocations = match.locations.findAll { loc ->
                loc.end - loc.start + 1 >= minResiduesLength // filter out locations with less than 8 residues
            }
            validLocations ? [(modelAccession): match] : null
        }.findAll { it != null }]
    }

    json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

def stockholmSeedParser(String pfamASeedFile) {
    Map<String, List<Map<String, List<String>>>> nestingInfo = [:]
    String accession = null
    new File(pfamASeedFile).eachLine { rawLine ->
        def line = decode(rawLine.bytes)
        if (line.startsWith("#=GF AC")) {
            accession = line.split("\\s+")[2].split("\\.")[0]
        } else if (line.startsWith("#=GF NE")) {
            String nestedAcc = line.split("\\s+")[2].replace(";", "")
            nestingInfo[accession]?.nested?.add(nestedAcc) ?:
                (nestingInfo[accession] = [nested: [nestedAcc]])
        }
    }
    return nestingInfo
}

def stockholmClansParser(String pfamClansFile, Map nestingInfo) {
    String accession = null
    new File(pfamClansFile).eachLine { rawLine ->
        def line = decode(rawLine.bytes)
        if (line.startsWith("#=GF AC")) {
            accession = line.split("\\s+")[2].split("\\.")[0]
        } else if (line.startsWith("#=GF MB")) {
            String modelAccession = line.split("\\s+")[2].replace(";", "")
            if (modelAccession) {
                if (!nestingInfo.containsKey(modelAccession)) {
                    nestingInfo[modelAccession] = [clan: accession]
                } else {
                    nestingInfo[modelAccession].clan = accession
                }
            }
        }
    }
    return nestingInfo
}

def stockholmDatParser(String pfamADatFile) {
    Map<String, String> name2acc = [:]
    Map<String, List<String>> acc2nested = [:]
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
            acc2nested[accession]?.add(nestedAcc) ?: (acc2nested[accession] = [nestedAcc] as Set)
        }
    }
    Map<String, List<String>> parsedDat = [:]
    acc2nested.each { acc, nestedNames ->
        nestedNames.each { nestedName ->
            String nestedAccession = name2acc[nestedName]
            parsedDat[acc]?.add(nestedAccession) ?: (parsedDat[acc] = [nestedAccession] as List)
        }
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

def createFragment(Match rawMatch, Map fragment, String fragmentDcStatus, boolean twoActualRegions) {
    int newLocationStart = rawMatch.locations[0].start
    int newLocationEnd = rawMatch.locations[0].end
    int finalLocationEnd = rawMatch.locations[0].end

    if (fragment['start'] <= newLocationStart && fragment['end'] >= newLocationEnd) {
        fragmentDcStatus = "NC_TERMINAL_DISC"
    } else if (fragmentDcStatus == "CONTINUOUS") {
        fragmentDcStatus = null
    }

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
    return [newLocationStart, newLocationEnd, finalLocationEnd, fragmentDcStatus, twoActualRegions]
}

def isOverlapping(location1, location2) {
    Math.max(location1.start, location2.start) <= Math.min(location1.end, location2.end)
}
