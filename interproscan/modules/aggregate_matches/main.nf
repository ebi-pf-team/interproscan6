import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.JsonEncoding
import com.fasterxml.jackson.core.JsonToken

process AGGREGATE_SEQS_MATCHES {
    label 'small'
    // Aggregates sequence meta data with the corresponding match data
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)
    val(nucleic)

    output:
    path("seq_matches_aggreg.json")

    exec:
    JsonFactory factory = new JsonFactory()
    ObjectMapper mapper = new ObjectMapper()
    JsonParser seqParser = factory.createParser(new File(seqsPath.toString()))
    JsonParser matchesParser = factory.createParser(new File(matchesPath.toString()))

    def seqMatchesAggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }

    assert seqParser.nextToken() == JsonToken.START_OBJECT
    def matchesInfo = mapper.readValue(matchesParser, Map)
    while (seqParser.nextToken() != JsonToken.END_OBJECT) {
        String seqId = seqParser.getCurrentName()
        seqParser.nextToken()
        def seqInfo = mapper.readValue(seqParser, Map)

        if (nucleic) {
            seqInfo.each { orf ->
                FastaSequence protSequence = FastaSequence.fromMap(orf)
                String protMD5 = protSequence.md5
                seqMatchesAggreg[protMD5].sequence = protSequence.sequence
                seqMatchesAggreg[protMD5].md5 = protMD5
                seqMatchesAggreg[protMD5].xref << ["name": orf.id + " " orf.description", "id": orf.id]
                if (seqMatchesAggreg[protMD5].translatedFrom == null) {
                    seqMatchesAggreg[protMD5].translatedFrom = []
                }
                seqMatchesAggreg[protMD5].translatedFrom << orf.translatedFrom // add nucleic seq metadata
                if (matchesInfo[orf.id]) {
                    seqMatchesAggreg[protMD5].matches = matchesInfo[orf.id]
                }
            }
        } else {
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
    mapper.writeValue(new File(outputFilePath.toString()), seqMatchesAggreg)
}

process AGGREGATE_ALL_MATCHES {
    label 'small'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    JsonFactory factory = new JsonFactory()
    ObjectMapper mapper = new ObjectMapper()
    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    JsonGenerator generator = factory.createGenerator(new File(outputFilePath.toString()), JsonEncoding.UTF8)
    generator.setCodec(mapper)
    generator.writeStartArray()

    seqMatches.each { file ->
        JsonParser parser = factory.createParser(new File(file.toString()))
        assert parser.nextToken() == JsonToken.START_OBJECT
        while (parser.nextToken() != JsonToken.END_OBJECT) {
            parser.nextToken()
            def info = mapper.readValue(parser, Map)
            generator.writeObject(info)
        }
        parser.close()
    }
    generator.writeEndArray()
    generator.close()
}
