import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process COMBINE_MATCHES {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path("combined_matches.json")

    exec:
    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- stops a seqMd5 appearing multiple times in the output
    members_matches.each { matchesPath ->
        if (matchesPath) {  // can be null if LOOKUP_MATCHES failed or was skipped
            new ObjectMapper().readValue(new File(matchesPath.toString()), Map).each { seqMd5, matches ->
                aggregatedMatches.computeIfAbsent(seqMd5, { [:] }).putAll(matches)
            }
        }
    }

    String outputFilePath = task.workDir.resolve("combined_matches.json")
    def json = JsonOutput.toJson(aggregatedMatches)
    new File(outputFilePath.toString()).write(json)
}