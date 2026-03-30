package uk.ac.ebi.interpro

import java.nio.file.Path
import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.databind.ObjectMapper


class ProcessCombine {
    static void run(List<Path> inputPaths, Path outputPath) {
        def mapper = new ObjectMapper()
        def aggregatedMatches = [:]  // seqMd5 -> (modelAcc -> match)

        inputPaths.each { path ->
            def data = path.newReader().withCloseable { reader -> 
                mapper.readValue(reader, Map)
            }
            data.each { seqMd5, matches ->
                aggregatedMatches.computeIfAbsent(seqMd5) { [:] }.putAll(matches)
            }
        }

        def jf = new JsonFactory()
        outputPath.withWriter { writer ->
            jf.createGenerator(writer).withCloseable { generator ->
                mapper.writeValue(generator, aggregatedMatches)
            }
        }
    }
}
