import groovy.json.JsonOutput
import groovy.json.JsonSlurper


process AGGREGATE_SEQS_MATCHES {
    // Aggregates sequence meta data with the corresponding match data
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)
    val(nucleic)

    output:
    path("seq_matches_aggreg.json")

    exec:
    JsonSlurper jsonSlurper = new JsonSlurper()
    def seqsInfo = jsonSlurper.parse(seqsPath)
    def matchesInfo = jsonSlurper.parse(matchesPath)

    def seqMatchesAggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }

    seqsInfo.each { seqKey, info ->
        /* seqKey = seq ID when input consists of protein seqs
           seqKey = md5 of ORF protein seq when input consists of nucleic seqs */

        String md5
        if (nucleic) {
            info.each { orf ->
                seqId = orf.id
                FastaSequence sequence = FastaSequence.fromMap(orf)
                md5 = sequence.md5
                seqMatchesAggreg[md5].sequence = sequence.sequence
                seqMatchesAggreg[md5].md5 = md5
                seqMatchesAggreg[md5].xref << ["name": "${sequence.id} ${sequence.description}", "id": sequence.id]
                seqMatchesAggreg[md5].translatedFrom = orf.translatedFrom // add nucleic seq metadata
            }
        } else {
            FastaSequence sequence = FastaSequence.fromMap(orf)
            md5 = sequence.md5
            seqId = seqKey
            seqMatchesAggreg[md5].sequence = sequence.sequence
            seqMatchesAggreg[md5].md5 = md5
            seqMatchesAggreg[md5].xref << ["name": "${sequence.id} ${sequence.description}", "id": sequence.id]
        }

        if (matchesInfo[seqId]) {  // the nucleic seq matches Map is keyed by the OrfId
            seqMatchesAggreg[md5].matches = matchesInfo[seqId]
        }
    }
    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    def json = JsonOutput.toJson(seqMatchesAggreg)
    new File(outputFilePath.toString()).write(json)
}

process AGGREGATE_ALL_MATCHES {
    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    def aggregatedResults = []
    JsonSlurper jsonSlurper = new JsonSlurper()
    seqMatches.each { file ->
        def jsonContent = jsonSlurper.parse(file)
        jsonContent.each { md5, info ->
            aggregatedResults << info
        }
    }

    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    def json = JsonOutput.toJson(aggregatedResults)
    new File(outputFilePath.toString()).write(json)
}
