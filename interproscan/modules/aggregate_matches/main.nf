import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ArrayNode
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.SerializationFeature

process AGGREGATE_SEQS_MATCHES {
    label 'local'
    // Aggregates sequence meta data with the corresponding match data
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)
    val(nucleic)

    output:
    path("seq_matches_aggreg.json")

    exec:
    /*
    Content of seqsPath when input is nucleic seqs:
      {prot seq md5: [{id: str, description: str, sequence: protein, md5: str, translatedFrom: {id, desc, seq (nucleic), md5, translatedFrom: null}]
    Content seqsPath when input is protein seqs:
      {seqId: {id: str, description: str, sequence: str, md5: str, translatedFrom: null}
    Content of matchesPath:
      {seqId: {sigAcc: {Match object represented as a Map}}
    Final output of this process:
      {prot seq md5: {Match object represented as a Map}}
    */
    // Build a single mapper for all readers and writers to save memory
    ObjectMapper jacksonMapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
    // Load the entire JSON file of matches to retrieve matches by the md5 of the seq
    def matchesMap = JsonReader.load(matchesPath.toString(), jacksonMapper)  // [seqId: matches]

    def seqMatchesAggreg = [:].withDefault { [sequence: '', md5: '', matches: [], xref: []] }

    JsonReader.stream(seqsPath.toString(), jacksonMapper) { JsonNode node ->
        if (nucleic) {  // node = [protMd5: [{protein}, {protein}, {protein}]]
            node.fields().each { entry ->
                JsonNode proteins = entry.value  // array of proteins from the predicted ORFs
                proteins.forEach { protein ->
                    processProteinData(protein, seqMatchesAggreg, matchesMap)
                    md5 = protein.get("md5").asText()
                    seqMatchesAggreg[md5].translatedFrom = seqMatchesAggreg[md5].translatedFrom ?: []
                    seqMatchesAggreg[md5].translatedFrom << protein.get("translatedFrom")
                }
            }
        } else {  // node = [protSeqId: {Seq}]
           processProteinData(node, seqMatchesAggreg, matchesMap)
        }
    }

    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    JsonWriter.writeMap(outputFilePath.toString(), jacksonMapper, seqMatchesAggreg)
}

def processProteinData(JsonNode protein, def seqMatchesAggreg, def matchesMap) {
    md5 = protein.get("md5").asText()
    seqId = protein.get("id").asText()
    seqMatchesAggreg[md5].sequence = protein.get("sequence").asText()
    seqMatchesAggreg[md5].md5 = md5
    seqMatchesAggreg[md5].xref << ["name": seqId + " " + protein.get("description"), "id": seqId]
    seqMatchesAggreg[md5].matches = matchesMap[seqId] ?: seqMatchesAggreg[md5].matches
}

process AGGREGATE_ALL_MATCHES {
    label 'local'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    // Build a single mapper for all readers and writers to save memory
    ObjectMapper jacksonMapper = new ObjectMapper().enabled(SerializationFeature.INDENT_OUTPUT)
    // Gather seqs by md5 to check for md5s with identical seqs and differing seq ids
    def allAggregatedData = [:]  // [md5: matches]
    seqMatches.each { file ->
        JsonReader.stream(file.toString(), jacksonMapper) { JsonNode node ->
            md5 = node.get("md5").asText()
            if (allAggregatedData.containsKey(md5)) {
                // merge existing and new matches
                ObjectNode existingNode = (ObjectNode) allAggregatedData.get(md5)          // mutable
                JsonNode newMatches = node.get("matches")                                  // immutable
                ((ObjectNode) existingNode.get("matches")).setAll((ObjectNode) newMatches)
                // merge existing and new xrefs - setAll is not applicable for ArrayNodes
                ArrayNode existingXref = (ArrayNode) existingNode.get("xref")
                ArrayNode newXref = (ArrayNode) node.get("xref")
                existingXref.addAll(newXref)
            } else {
                allAggregatedData[md5] = node
            }
        }
    }

    // write the output as an array of seq objects [{Seq}, {Seq}, {Seq}]
    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    JsonWriter.stream("${file_path}", mapper) { generator ->
        allAggregatedData.each { md5, jsonNode ->
            generator.writeTree(jsonNode)
        }
    }
}
