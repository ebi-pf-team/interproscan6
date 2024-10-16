import groovy.json.JsonOutput

process RUN_SIGNALP {
    label 'signalp_runner'

    input:
    tuple val(meta), path(fasta)
    val orgType
    val mode

    output:
    tuple val(meta), path("signalp_out"), val(orgType)

    script:
    """
    signalp6 \
        --organism ${orgType} \
        --fastafile ${fasta} \
        --output_dir signalp_out \
        --mode ${mode}
    """
}

process PARSE_SIGNALP {
    label 'analysis_parser'

    input:
    tuple val(meta), val(signalp_out), val(orgType)
    val threshold

    output:
    tuple val(meta), path("signalp.json")

    exec:
    /*
    Don't parse only the JSON as we want the start and end positions
    that are annotated by SignalP for the rare cases were the signal
    peptide does not start at position 1.
    */
    String signalDir = signalp_out.toString()
    String gff3FileName = "output.gff3"
    String txtFileName = "prediction_results.txt"
    float threshold = threshold as float
    String orgType = orgType as String
    def Map<String, Match> hits = new LinkedHashMap<>()

    // Retrieve all signal peptide hits from the gff3 file
    File gff3File = new File(signalDir, gff3FileName)
    gff3File.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null) {
            def matcher = line =~ ~/^(.+?)\sSignalP-\d+\.\d+\ssignal_peptide\s+(\d+)\s(\d+)\s(\d+\.?\d*)\s\.\s\.\s\.$/
            if (matcher.find()) {
                String queryAccession = matcher.group(1).trim()
                float pvalue = Float.parseFloat(matcher.group(4).trim())
                if (pvalue >= threshold) {
                    Match match = new Match("SignalP")
                    Location location = new Location(
                        Integer.parseInt(matcher.group(2).trim()),  // start
                        Integer.parseInt(matcher.group(3).trim()),  // end
                        pvalue
                    )
                    match.addLocation(location)
                    hits.put(queryAccession, match)
                }
            }
        }
    }

    // Retrieve the cleavage site start/end from the TSV (.txt) file, and add to the matches
    File tsvFile = new File(signalDir, txtFileName)
    tsvFile.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null) {
            def matcher = line =~ ~/^(.+?)\sSP\s.+?CS pos: (\d+)-(\d+)\.\sPr:\s\d+\.?\d+$/
            if (matcher.find()) {
                String queryAccession = matcher.group(1).trim()
                if (hits.containsKey(queryAccession)) {  // protein may not have passed the earlier threshold check
                    hits[queryAccession].addSignalPeptide(
                            orgType,
                            matcher.group(2).trim().toInteger(),  // cleavageSite Start
                            matcher.group(3).trim().toInteger()   // cleavageSite End
                    )
                }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("signalp.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
