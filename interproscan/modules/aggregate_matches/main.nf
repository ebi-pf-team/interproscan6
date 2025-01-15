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
    def outputFilePath = task.workDir.resolve("seq_matches_aggreg.json")
    JsonGenerator generator = factory.createGenerator(new File(outputFilePath.toString()), JsonEncoding.UTF8)
    generator.setCodec(mapper)
    generator.writeStartObject()
    JsonParser seqParser = factory.createParser(new File(seqsPath.toString()))
    JsonParser matchesParser = factory.createParser(new File(matchesPath.toString()))
    assert seqParser.nextToken() == JsonToken.START_OBJECT
    def matchesMap = mapper.readValue(matchesParser, Map)

    while (seqParser.nextToken() != JsonToken.END_OBJECT) {
        String seqId = seqParser.getCurrentName()
        seqParser.nextToken()
        def seqInfo = mapper.readValue(seqParser, Map)
        if (nucleic) {
            seqInfo.each { orf ->
                FastaSequence protSequence = FastaSequence.fromMap(orf)
                String protMD5 = protSequence.md5

                generator.writeFieldName(protMD5)
                generator.writeStartObject()
                generator.writeStringField("sequence", protSequence.sequence)
                generator.writeStringField("md5", protMD5)
                generator.writeArrayFieldStart("xref")
                generator.writeStartObject()
                generator.writeStringField("name", "${orf.id} ${orf.description}")
                generator.writeStringField("id", orf.id)
                generator.writeEndObject()
                generator.writeEndArray()

                if (orf.translatedFrom != null) {
                    generator.writeArrayFieldStart("translatedFrom")
                    generator.writeString(orf.translatedFrom)
                    generator.writeEndArray()
                }

                def matchData = matchesMap[orf.id] ?: []
                generator.writeObjectFieldStart("matches")
                matchData.each { key, value ->
                    generator.writeFieldName(key)
                    generator.writeObject(value)
                }
                generator.writeEndObject()
            }
        } else {
            FastaSequence sequence = FastaSequence.fromMap(seqInfo)
            String md5 = sequence.md5

            generator.writeFieldName(md5)
            generator.writeStartObject()
            generator.writeStringField("sequence", sequence.sequence)
            generator.writeStringField("md5", md5)
            generator.writeArrayFieldStart("xref")
            generator.writeStartObject()
            generator.writeStringField("name", "${sequence.id} ${sequence.description}")
            generator.writeStringField("id", sequence.id)
            generator.writeEndObject()
            generator.writeEndArray()

            def matchData = matchesMap[seqId] ?: []
            generator.writeObjectFieldStart("matches")
            matchData.each { key, value ->
                generator.writeFieldName(key)
                generator.writeObject(value)
            }
            generator.writeEndObject()
        }
    }
    seqParser.close()
    generator.writeEndObject()
    generator.close()
}

process AGGREGATE_ALL_MATCHES {
    label 'small'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    ObjectMapper mapper = new ObjectMapper()
    JsonFactory factory = new JsonFactory()
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

