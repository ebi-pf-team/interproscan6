import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process COMBINE_MATCHES {
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path("combined_matches.json")

    exec:
    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- stops a seqMd5 appearing multiple times in the output
    members_matches.each{ matchesPath ->
        def matchesFileMap = new ObjectMapper().readValue(new File(matchesPath.toString()), Map)
        matchesFileMap.each { String seqMd5, Map matches ->
            def seqEntry = aggregatedMatches.computeIfAbsent(seqMd5, { [:] })
            matches.each { String modelAcc, Map match ->
                seqEntry[modelAcc] = match
            }
        }
    }

    String outputFilePath = task.workDir.resolve("combined_matches.json")
    def json = JsonOutput.toJson(aggregatedMatches)
    new File(outputFilePath.toString()).write(json)
}