import groovy.json.JsonSlurper
import com.fasterxml.jackson.databind.ObjectMapper


process AGGREGATE_SEQS_MATCHES {
    label 'small'
    // Aggregates sequence meta data with the corresponding match data
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)
    val(nucleic)

    output:
    path("seq_matches_aggreg.json")

    exec:
    JsonSlurper jsonSlurper = new JsonSlurper()
    def mapper = new ObjectMapper()
    def seqsInfo = jsonSlurper.parse(seqsPath)
    def matchesInfo = jsonSlurper.parse(matchesPath)

    def seqMatchesAggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }

    seqsInfo.each { seqId, info ->
        if (nucleic) {
            info.each { orf ->
                FastaSequence protSequence = FastaSequence.fromMap(orf)
                protMD5 = protSequence.md5
                seqMatchesAggreg[protMD5].sequence = protSequence.sequence
                seqMatchesAggreg[protMD5].md5 = protMD5
                seqMatchesAggreg[protMD5].xref << ["name": "${orf.id} ${orf.description}", "id": orf.id]
                if (seqMatchesAggreg[protMD5].translatedFrom == null) {
                    seqMatchesAggreg[protMD5].translatedFrom = []
                }
                seqMatchesAggreg[protMD5].translatedFrom << orf.translatedFrom // add nucleic seq metadata
                if (matchesInfo[orf.id]) {
                    seqMatchesAggreg[protMD5].matches = matchesInfo[orf.id]
                }
            }
        } else {
            FastaSequence sequence = FastaSequence.fromMap(info)
            md5 = sequence.md5
            seqMatchesAggreg[md5].sequence = sequence.sequence
            seqMatchesAggreg[md5].md5 = md5
            seqMatchesAggreg[md5].xref << ["name": "${sequence.id} ${sequence.description}", "id": sequence.id]
            if (matchesInfo[seqId]) {
                seqMatchesAggreg[md5].matches = matchesInfo[seqId]
            }
        }
    }
    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    mapper.writeValue(new File(outputFilePath.toString()), seqMatchesAggreg)    
}

process AGGREGATE_ALL_MATCHES {
    label 'small'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    def aggregatedResults = []
    JsonSlurper jsonSlurper = new JsonSlurper()
    def mapper = new ObjectMapper()

    seqMatches.each { file ->
        def jsonContent = jsonSlurper.parse(file)
        jsonContent.each { md5, info ->
            aggregatedResults << info
        }
    }

    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    mapper.writeValue(new File(outputFilePath.toString()), aggregatedResults)    
}
