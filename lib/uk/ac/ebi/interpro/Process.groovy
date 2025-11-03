package uk.ac.ebi.interpro

import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

class Process {
    static void combineMatches(List<String> inputPaths, String outputPath) {
        def mapper = new ObjectMapper()
        def aggregatedMatches = [:]  // seqMd5 -> (modelAcc -> match)

        inputPaths.each { path ->
            def file = new File(path)
            def data = mapper.readValue(file, Map)
            data.each { seqMd5, matches ->
                aggregatedMatches.computeIfAbsent(seqMd5) { [:] }.putAll(matches)
            }
        }

        def json = JsonOutput.toJson(aggregatedMatches)
        new File(outputPath).text = json
    }
}
