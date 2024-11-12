import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process AGGREGATE_SEQS_MATCHES {
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)

    output:
    path("seq_matches_aggreg.json")

    exec:
    JsonSlurper jsonSlurper = new JsonSlurper()
    def seqsInfo = jsonSlurper.parse(seqsPath)
    def matchesInfo = jsonSlurper.parse(matchesPath)

    def seq_matches_aggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }
    seqsInfo.each { seqId, info ->
        FastaSequence sequence = FastaSequence.fromMap(info)
        String md5 = sequence.md5
        seq_matches_aggreg[md5].sequence = sequence.sequence
        seq_matches_aggreg[md5].md5 = md5
        seq_matches_aggreg[md5].xref << ["name": "${sequence.id} ${sequence.description}", "id": sequence.id]
        if (matchesInfo[seqId]) {
            seq_matches_aggreg[md5].matches = matchesInfo[seqId]
        }
    }
    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    def json = JsonOutput.toJson(seq_matches_aggreg)
    new File(outputFilePath.toString()).write(json)
}

process AGGREGATE_RESULTS {
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
