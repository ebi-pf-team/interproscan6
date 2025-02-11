import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.SerializationFeature

process AGGREGATE_MATCHES {
    // For each sequence, we aggregate the matches from all applications
    label 'local'

    input:
    val matchesFiles

    output:
    path("aggregated_results.json")

    exec:
    // Build a single mapper for all readers and writers to save memory
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)

    // Gather seqs by md5 to check for md5s with identical seqs and differing seq ids
    def aggregatedMatches = [:]
    matchesFiles.each { file ->    // seqMatches = [md5: Map<seq: str, md5: str, matches: {obj}|[], xref: [{obj}]>]]
        JsonReader.streamJson(file.toString(), jacksonMapper) { String md5, JsonNode node ->
            if (aggregatedMatches.containsKey(md5)) {
                // merge existing and new matches
                ObjectNode existingMatches = (ObjectNode) aggregatedMatches.get("matches")  // mutable
                JsonNode newMatches = node.get("matches")                                   // immutable
                ((ObjectNode) existingMatches.get("matches")).setAll((ObjectNode) newMatches)
            } else {
                aggregatedMatches[md5] = node
            }
        }
    }

    // write the output as an array of seq objects [{Seq}, {Seq}, {Seq}]
    def outputFilePath = task.workDir.resolve("aggregated_matches.json")
    JsonWriter.streamArray(outputFilePath.toString(), jacksonMapper) { generator ->
        aggregatedMatches.each { md5, jsonNode ->
            JsonWriter.writeMap(generator, jacksonMapper, jsonNode)
        }
    }
}
