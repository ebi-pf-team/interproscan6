import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import java.util.regex.Pattern

process RUN_PFAM {
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
    val minLength
    val seedPath
    val clanPath
    val datPath

    output:
    tuple val(meta), path("pfam.json")

    exec:
    def outputFilePath = task.workDir.resolve("pfam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString())

    nestedSeeds = stockholmSeedParser(seedPath)
    clans = stockholmClansParser(clanPath, nestedSeeds)
    dat = stockholmDatParser(datPath)

    filteredMatches = [:]
    hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
        // sorting matches by evalue ASC, score DESC
        def sortedMatches = matches.sort { a, b ->
            def evalueComparison = a.value.evalue <=> b.value.evalue
            def scoreComparison = -(a.value.score <=> b.value.score)
            evalueComparison != 0 ? evalueComparison : scoreComparison
        }

        filteredMatches[seqId] = [:]
        sortedMatches.each { modelAccession, match ->
            boolean keep = true
            modelAccession = modelAccession.split("\\.")[0]
            def candidateMatch = clans[modelAccession] ?: [:]
            def candidateClan = candidateMatch?.get("clan", null)
            filteredMatches[seqId].each { filteredAcc, filteredMatch ->
                filteredMatch = clans[filteredAcc] ?: [:]
                def filteredClan = filteredMatch?.get("clan", null)
                if (candidateClan && candidateClan == filteredClan) {  // enter if same clan
                    def notOverlappingLocations = []
                    match.locations.each { candidateLocation ->
                        boolean overlapping = false
                        filteredMatches.locations.each { filteredLocation ->
                            if (Math.max(candidateLocation.start, filteredLocation.start) <= Math.min(candidateLocation.end, filteredLocation.end)) {  // checking overlapping locations
                                def candidateNested = candidateMatch?.get("nested", [])
                                def filteredNested = filteredMatch?.get("nested", [])
                                matchesAreNested = (candidateNested != null && filteredNested != null) && candidateNested.intersect(filteredNested).size() > 0
                                if (!matchesAreNested) {  // case of overlapped but NOT nested
                                    overlapping = true
                                    return
                                }
                            }
                        }
                        if (!overlapping) notOverlappingLocations << candidateLocation
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

    // building fragments
    processedMatches = [:]
    filteredMatches.collectEntries { seqId, matches ->
        matches.each { modelAccession, match ->
            if (!match.locations) {
                return
            }
            def nestedModels = dat.get(modelAccession, [])
            if (nestedModels) {
                List locationFragments = []
                matches.each { otherAccession, otherMatch ->
                    otherMatch.locations.each { location ->
                        overlapping = Math.max(location.start, match.locations[0].start) <= Math.min(location.end, match.locations[0].end)
                        if (otherAccession in nestedModels && overlapping) {
                            locationFragments << [
                                start: location.start,
                                end: location.end
                            ]
                        }
                    }
                }
                locationFragments.sort { a, b ->
                    a['start'] <=> b['start'] ?: a['end'] <=> b['end']
                }
                def fragmentDcStatus = "CONTINUOUS"
                List rawDiscontinuousMatches = [match]
                locationFragments.each { fragment ->
                    boolean twoActualRegions = false
                    (newLocationStart, newLocationEnd, finalLocationEnd,
                        fragmentDcStatus, twoActualRegions) = createFragment(
                            rawDiscontinuousMatches[0],
                            fragment,
                            fragmentDcStatus,
                            twoActualRegions
                    )
                    rawDiscontinuousMatches[0].locations.each { loc ->
                        loc.fragments = []
                        loc.fragments << new LocationFragment(
                            newLocationStart,
                            newLocationEnd,
                            fragmentDcStatus)
                        if (twoActualRegions) {
                            loc.fragments << new LocationFragment(
                                fragment['end'] + 1,
                                finalLocationEnd,
                                "N_TERMINAL_DISC")
                        }
                    }
                }

                rawDiscontinuousMatches.each { rawDiscontinuousMatch ->
                    rawDiscontinuousMatch.locations.each { location ->
                        def matchLength = location.end - location.start + 1
                        if (matchLength >= minLength) {
                            processedMatches.putIfAbsent(seqId, [:])
                            processedMatches[seqId].putIfAbsent(modelAccession, [])
                            processedMatches[seqId][modelAccession] = rawDiscontinuousMatch
                        }
                    }
                }
            } else {
                if (match.locations.size() > 0) {
                    def matchLength = (match.locations[0].end) - (match.locations[0].start) + 1
                    if (matchLength >= minLength) {
                        processedMatches.putIfAbsent(seqId, [:])
                        processedMatches[seqId][modelAccession] = match
                    }
                }
            }
        }
    }

    json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

def stockholmSeedParser(String pfamASeedFile) {
    Map nestingInfo = [:]
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
                nestingInfo[modelAccession]?.clan = accession ?:
                    (nestingInfo[modelAccession] = [clan: accession])
            }
        }
    }
    return nestingInfo
}

def stockholmDatParser(String pfamADatFile) {
    Map name2acc = [:]
    Map acc2nested = [:]
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
    Map parsedDat = [:]
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

    if ((fragment['start']) <= newLocationStart && (fragment['end']) >= newLocationEnd) {
        fragmentDcStatus = "NC_TERMINAL_DISC"
    } else if (fragmentDcStatus == "CONTINUOUS") {
        fragmentDcStatus = null
    }

    if ((fragment['start']) <= newLocationStart) {
        newLocationStart = (fragment['end']) + 1
        fragmentDcStatus = "N_TERMINAL_DISC"
    } else if ((fragment['end']) >= newLocationEnd) {
        newLocationEnd = (fragment['start']) - 1
        fragmentDcStatus = "C_TERMINAL_DISC"
    } else if ((fragment['start']) > newLocationStart && (fragment['end']) < newLocationEnd) {
        newLocationEnd = (fragment['start']) - 1
        twoActualRegions = true
        fragmentDcStatus = "C_TERMINAL_DISC"
    }

    return [newLocationStart, newLocationEnd, finalLocationEnd, fragmentDcStatus, twoActualRegions]
}
