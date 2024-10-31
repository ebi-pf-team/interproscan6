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

    // Return Pfam clans AND nesting relationships between models
    seedNesting = stockholmSeedParser(seedPath)
    clan = stockholmClansParser(clanPath, seedNesting)
    dat = stockholmDatParser(datPath)

//     hmmerMatches = hmmerMatches.collectEntries { seqId, matches ->
//         matches.each { modelAccession, match ->
//
//         }
//         return [(seqId): matches]
//     }
//
//     def json = JsonOutput.toJson(hmmerMatches)
//     new File(outputFilePath.toString()).write(json)
}

def stockholmSeedParser(String pfamASeedFile) {
    def nestingInfo = [:]
    def accession = null
    new File(pfamASeedFile).withInputStream { stream ->
        stream.eachLine { rawLine ->
            def line = decode(rawLine.bytes)
            if (line.startsWith("#=GF AC")) {
                accession = line.split("\\s+")[2]
            } else if (line.startsWith("#=GF NE")) {
                def nestedAcc = line.split("\\s+")[2].replace(";", "")
                nestingInfo[accession]?.nestedAcc?.add(nestedAcc) ?:
                    (nestingInfo[accession] = [nested: [nestedAcc]])
            }
        }
    }
//     println "nestingInfo: $nestingInfo"
    return nestingInfo
}

def stockholmClansParser(String pfamClansFile, Map parsedSeed) {
    def parsedResult = [:]
    def accession = null
    new File(pfamClansFile).withInputStream { stream ->
        stream.eachLine { rawLine ->
            def line = decode(rawLine.bytes)
            if (line.startsWith("#=GF AC")) {
                accession = line.split("\\s+")[2]
            } else if (line.startsWith("#=GF MB")) {
                def modelAccession = line.split("\\s+")[2].replace(";", "")
                if (modelAccession) {
                    parsedResult[modelAccession]?.clan = accession ?:
                        (parsedResult[modelAccession] = [clan: accession])
                }
            }
        }
    }
    println "Parsed clans: $parsedResult"
    return parsedResult
}

def stockholmDatParser(String pfamADatFile) {
    def name2acc = [:]
    def acc2nested = [:]
    name = null
    accession = null
    new File(pfamADatFile).withInputStream { stream ->
        stream.eachLine { rawLine ->
            def line = decode(rawLine.bytes)
            if (line.startsWith("#=GF ID")) {
                def name = line.split()[2]
            } else if (line.startsWith("#=GF AC")) {
                def accession = line.split()[2].split("\\.")[0]
                name2acc[name] = accession
            } else if (line.startsWith("#=GF NE")) {
                def nestedAcc = line.split()[2]
                acc2nested[accession]?.add(nestedAcc) ?: (acc2nested[accession] = [nestedAcc] as Set)
            }
        }
    }
    def parsedDat = [:]
    acc2nested.each { acc, nestedNames ->
        nestedNames.each { name ->
            def nestedAccession = name2acc[name]
            parsedDat[acc]?.add(nestedAccession) ?: (parsedDat[acc] = [nestedAccession] as List)
        }
    }
    println "Parsed dat: $parsedDat"
    return parsedDat
}

def decode(byte[] b) {
    try {
        return new String(b, "UTF-8").trim()
    } catch (Exception e) {
        return new String(b, "ISO-8859-1").trim()
    }
}
