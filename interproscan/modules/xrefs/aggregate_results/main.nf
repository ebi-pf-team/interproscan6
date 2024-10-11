import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process AGGREGATE_RESULTS {
    label 'xref'

    input:
    tuple val(meta), val(matches2xrefs)

    output:
    path "aggregated_matches2xrefs.json"

    exec:
    def aggregatedXrefs = [:]
    JsonSlurper jsonSlurper = new JsonSlurper()
    matches2xrefs.each { filePath ->
        def fileName = filePath.tokenize('/').last()
        def matches = jsonSlurper.parse(new File(filePath)).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                Match matchObject = Match.fromMap(jsonMatch)
                if (fileName.contains("matches2entries")) {
                    match = matchObject
                } else if (fileName.contains("matches2go")) {
                    goTerms = matchObject.signature.entry.goXrefs ?: []
                    match.signature.entry.goXrefs = goTerms
                } else if (fileName.contains("matches2pa")) {
                    pathways = matchObject.signature.entry.pathwaysXrefs ?: []
                    match.signature.entry.pathwaysXrefs = pathways
                }
                aggregatedXrefs[seqId] = match
            }]
        }
    }
    def outputFilePath = task.workDir.resolve("aggregated_matches2xrefs.json")
    def json = JsonOutput.toJson(aggregatedXrefs)
    new File(outputFilePath.toString()).write(json)
}
