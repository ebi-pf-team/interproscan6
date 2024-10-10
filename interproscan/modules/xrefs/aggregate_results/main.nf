import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process AGGREGATE_RESULTS {
    label 'xref'

    input:
    tuple val(meta), val(matches2entries)
    tuple val(meta), val(xrefsMix)

    output:
    path "aggregated_matches2xrefs.json"

    exec:
    def matchObjectAggregated = [:]
    JsonSlurper jsonSlurper = new JsonSlurper()
    matches = jsonSlurper.parse(xrefsMix).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            Match matchObject = Match.fromMap(jsonMatch)
            if (!matchObjectAggregated.containsKey(seqId)) {
                matchObjectAggregated[seqId] = []
            }
            matchObjectAggregated[seqId].add(matchObject)
            [(matchId): jsonMatch]
        }]
    }
    def outputFilePath = task.workDir.resolve("aggregated_matches2xrefs.json")
    def json = JsonOutput.toJson(matchObjectAggregated)
    new File(outputFilePath.toString()).write(json)
}
