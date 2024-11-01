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
    val min_length
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
            accession = modelAccession.split("\\.")[0]
            def candidateMatch = clans[accession] ?: [:]
            def candidateClan = candidateMatch?.get("clan", null)
            filteredMatches[seqId].each { filteredAcc, filteredMatch ->
                filteredMatch = clans[filteredAcc] ?: [:]
                def filteredClan = filteredMatch?.get("clan", null)
                if (candidateClan && candidateClan == filteredClan) {  // case of same clan
                    println "Same clan: ${candidateClan}"
                    def notOverlappingLocations = []
                    match.locations.each { candidateLocation ->  // checking overlapping locations
                        boolean overlapping = false
                        filteredMatches.locations.each { filteredLocation ->
                            if (matchesOverlap(candidateLocation, filteredLocation)) {
                                def candidateNested = candidateMatch?.get("nested", [])
                                def filteredNested = filteredMatch?.get("nested", [])
                                if (!matchesAreNested(candidateNested, filteredNested)) {  // case of overlapped but NOT nested
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

    finalMatches = filteredMatches.collectEntries { seqId, matches ->
        matches.each { modelAccession, match ->
            if (!match.locations) {
                return
            }
            String accession = modelAccession.split("\\.")[0]
            def nestedModels = dat.get(accession, [])
            //build fragments...
        }
    }

    json = JsonOutput.toJson(finalMatches)
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

def matchesOverlap(Map one, Map two) {
    return Math.max(one['start'].toInteger(), two['start'].toInteger()) <= Math.min(one['end'].toInteger(), two['end'].toInteger())
}

def matchesAreNested(List one, List two) {
    if (one == null || two == null) {
        return false
    }
    return one.intersect(two).size() > 0
}
