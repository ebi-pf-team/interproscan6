package uk.ac.ebi.interpro

import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper

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

        def jf = new JsonFactory()
        new File(outputPath).withWriter { writer ->
            def gen = jf.createGenerator(writer)
            mapper.writeValue(gen, aggregatedMatches)
            gen.close()
        }
    }
}
