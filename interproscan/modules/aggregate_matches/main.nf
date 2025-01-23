import com.fasterxml.jackson.core.JsonToken

process AGGREGATE_SEQS_MATCHES {
    label 'local'
    // Aggregates sequence meta data with the corresponding match data
    input:
    tuple val(meta), val(seqsPath), val(matchesPath)
    val(nucleic)

    output:
    path("seq_matches_aggreg.json")

    exec:
    JsonProcessor processor = new JsonProcessor()
    def seqParser = processor.createParser(seqsPath.toString())
    def matchesParser = processor.createParser(matchesPath.toString())

    def seqMatchesAggreg = [:].withDefault { [
        sequence: '',
        md5: '',
        matches: [],
        xref: []
    ] }

    def matchesInfo = processor.jsonToMap(matchesParser)
    while (seqParser.nextToken() != JsonToken.END_OBJECT) {
        String seqId = seqParser.getCurrentName()
        seqParser.nextToken()
        if (nucleic) {
            def seqInfo = processor.jsonToList(seqParser)
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
            def seqInfo = processor.jsonToMap(seqParser)
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
    processor.write(outputFilePath.toString(), seqMatchesAggreg)
}

process AGGREGATE_ALL_MATCHES {
    label 'local'

    input:
    val seqMatches

    output:
    path("aggregated_results.json")

    exec:
    JsonProcessor processor = new JsonProcessor()
    def outputFilePath = task.workDir.resolve("aggregated_results.json")
    def generator = processor.createGenerator(outputFilePath.toString())

    generator.writeStartArray()
    seqMatches.each { file ->
        def parser = processor.createParser(file.toString())
        while (parser.nextToken() != JsonToken.END_OBJECT) {
            parser.nextToken()
            def info = processor.jsonToMap(parser)
            generator.writeObject(info)
        }
        parser.close()
    }
    generator.writeEndArray()
    generator.close()
}
