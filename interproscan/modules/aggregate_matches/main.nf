import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.JsonNode
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
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    // Load the entire JSON file of matches to retrieve matches by the md5 of the seq
    def matchesMap = JsonReader.load(matchesPath.toString(), jacksonMapper)  // [seqId: [modelAcc: [Match]]

    def seqMatchesAggreg = [:]
    JsonReader.streamJson(seqsPath.toString(), jacksonMapper) { String seqId, JsonNode node ->
        if (nucleic) { // seqId = Protein Seq MD5, node = [{prot}, {prot}]
            node.forEach { protein ->
                protSeqId = protein.get("id").asText()
                processProteinData(protein, seqMatchesAggreg, matchesMap, protSeqId)
                seqMatchesAggreg[seqId].translatedFrom = seqMatchesAggreg[seqId].translatedFrom ?: []
                seqMatchesAggreg[seqId].translatedFrom << protein.get("translatedFrom")
            }
        } else {  // node = [Protein Seq Id: {protein}]
            processProteinData(node, seqMatchesAggreg, matchesMap, seqId)
        }
    }

    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    JsonWriter.writeMaptoFile(outputFilePath.toString(), jacksonMapper, seqMatchesAggreg)
}

def processProteinData(JsonNode protein, Map seqMatchesAggreg,  Map<String, JsonNode> matchesMap, String seqId) {
    md5 = protein.get("md5").asText()
    def entry = seqMatchesAggreg.computeIfAbsent(md5, { [sequence: protein.get("sequence").asText(), md5: md5, matches: [:], xref: []] })
    entry.xref << ["name": seqId + " " + protein.get("description"), "id": seqId]
    if (matchesMap.containsKey(seqId)) { // matchesMap[seqId] is an ObjectNode, keyed by modelAcc and valued by ObjectNode repr of Matches
        matchesMap[seqId].fields().each { matchNode ->
            seqMatchesAggreg[md5]["matches"][matchNode.key] = Match.fromJsonNode((JsonNode) matchNode.value)
        }
    }
}

process AGGREGATE_ALL_MATCHES {
    label 'local'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    // Build a single mapper for all readers and writers to save memory
    ObjectMapper jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    // Gather seqs by md5 to check for md5s with identical seqs and differing seq ids
    def allAggregatedData = [:]
    seqMatches.each { file ->    // seqMatches = [md5: Map<seq: str, md5: str, matches: {obj}|[], xref: [{obj}]>]]
        JsonReader.streamJson(file.toString(), jacksonMapper) { String md5, JsonNode node ->
            if (allAggregatedData.containsKey(md5)) {
                // merge existing and new matches
                ObjectNode existingMatches = (ObjectNode) allAggregatedData.get("matches")  // mutable
                JsonNode newMatches = node.get("matches")                                   // immutable
                ((ObjectNode) existingMatches.get("matches")).setAll((ObjectNode) newMatches)
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
    JsonWriter.streamArray(outputFilePath.toString(), jacksonMapper) { generator ->
        allAggregatedData.each { md5, jsonNode ->
            JsonWriter.writeMap(generator, jacksonMapper, jsonNode)
        }
    }
}
