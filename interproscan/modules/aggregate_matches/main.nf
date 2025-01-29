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
    JsonProcessor jsonProcessor =
    def seqParser = JsonProcessor.createParser(seqsPath.toString())
    def matchesParser = JsonProcessor.createParser(matchesPath.toString())

    def seqMatchesAggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }

    def matchesInfo = JsonProcessor.jsonToMap(matchesParser)
    while (seqParser.nextToken() != JsonToken.END_OBJECT) {
        String seqId = seqParser.getCurrentName()
        seqParser.nextToken()
        if (nucleic) {
            def seqInfo = JsonProcessor.jsonToList(seqParser)
            seqInfo.each { orf ->
                FastaSequence protSequence = FastaSequence.fromMap(orf)
                String protMD5 = protSequence.md5
                seqMatchesAggreg[protMD5].sequence = protSequence.sequence
                seqMatchesAggreg[protMD5].md5 = protMD5
                seqMatchesAggreg[protMD5].xref << ["name": orf.id + " " + orf.description, "id": orf.id]
                if (seqMatchesAggreg[protMD5].translatedFrom == null) {
                    seqMatchesAggreg[protMD5].translatedFrom = []
                }
                seqMatchesAggreg[protMD5].translatedFrom << orf.translatedFrom // add nucleic seq metadata
                if (matchesInfo[orf.id]) {
                    seqMatchesAggreg[protMD5].matches = matchesInfo[orf.id]
                }
            }
        } else {
            def seqInfo = JsonProcessor.jsonToMap(seqParser)
            FastaSequence sequence = FastaSequence.fromMap(seqInfo)
            String md5 = sequence.md5
            seqMatchesAggreg[md5].sequence = sequence.sequence
            seqMatchesAggreg[md5].md5 = md5
            seqMatchesAggreg[md5].xref << ["name": sequence.id + " " + sequence.description, "id": sequence.id]
            if (matchesInfo[seqId]) {
                seqMatchesAggreg[md5].matches = matchesInfo[seqId]
            }
        }
    }

    seqParser.close()
    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    JsonProcessor.write(outputFilePath.toString(), seqMatchesAggreg)
}

process AGGREGATE_ALL_MATCHES {
    label 'local'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
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

    // write the output [{Seq}, {Seq}, {Seq}]
    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    JsonWriter.stream("${file_path}", mapper) { generator ->
        allAggregatedData.each { md5, jsonNode ->
            generator.writeTree(jsonNode)
        }
    }
}
