import groovy.json.JsonOutput

process RUN_PIRSF {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out"), path("hmmsearch.dtbl")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        -o hmmsearch.out \
        --domtblout hmmsearch.dtbl \
        ${hmmdb} ${fasta}
    """
}

process PARSE_PIRSF {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmsearch_out), val(hmmsearch_dtbl)
    val pirsf_dat_file

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def datData = PirsfDatEntry.parsePirsfDatFile(pirsf_dat_file.toString())  // [datEntries, datChildren]
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString())

    /* Retrieve the matches from the dtbl file because we need the tlen
    for the hmmLength, and PIRSF does not assign the same same hmmer start/ends
    to the same Match.Location fields as other member database. For example,
    PIRSF uses the envelope start and end for the location start/end. */
    Map<String, Map<String, Match>> matches = new LinkedHashMap<>()
    File hmmerDtblFile = new File(hmmsearch_dtbl.toString())
    hmmerDtblFile.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("#")) {
                def lineData = line.trim().split("\\s+")
                String proteinAccession = lineData[0]
                String modelAccession = lineData[4]

                if (!matches.containsKey(proteinAccession)) {
                    matches.put(proteinAccession, new LinkedHashMap<>())
                }

                if (!matches[proteinAccession].containsKey(modelAccession)) {
                    Match match = new Match(
                        modelAccession,
                        lineData[6].toDouble(),  // seqEvalue
                        lineData[7].toDouble(),  // seqScore
                        lineData[8].toDouble()   // seqBias
                    )
                    match.sequenceLength = lineData[5].toInteger()
                    matches[proteinAccession].put(modelAccession, match)
                }

                Location location = new Location(
                    lineData[19].toInteger(),   // Start: use the envelope
                    lineData[20].toInteger(),   // End: use the envelope
                    lineData[17].toInteger(),   // hmmStart: use ali from
                    lineData[18].toInteger(),   // hmmEnd: use ali to
                    lineData[2].toInteger(),    // hmmLength
                    null,                       // hmmBounds -- added below
                    lineData[19].toInteger(),   // envelopeStart
                    lineData[20].toInteger(),   // envelopeEnd
                    lineData[12].toDouble(),    // domain e-value
                    lineData[13].toFloat(),     // domain score
                    lineData[14].toFloat()      // domain bias
                )

                for (Match hmmerMatch: hmmerMatches[proteinAccession][modelAccession]) {
                    for (Location hmmerLocation in hmmerMatch.locations) {
                        if (location.start == hmmerLocation.envelopeStart &&
                            location.end == hmmerLocation.envelopeEnd) {
                            location.hmmBounds = hmmerLocation.hmmBounds
                        }
                    }
                }
                matches[proteinAccession][modelAccession].addLocation(location)
            }
        }
    }

    Map<String, Map<String, Match>> processedMatches = new LinkedHashMap<>()
    def datEntries = datData[0]
    def datChildren = datData[1]

    double LENGTH_RATIO_THRESHOLD = 0.67
    double OVERLAP_THRESHOLD = 0.8
    double LENGTH_DEVIATION_THRESHOLD = 3.5
    int MINIMUM_LENGTH_DEVIATION = 50

    def promote = [:]
    def store = [:]

    /* Filter the matches using the data in the
    dat file */
    matches.each { proteinAccession, modelMatches ->
        modelMatches.each { modelAccession, match ->
            match.locations.each { location ->
                // Calculate ratios
                double ovl = Math.abs(location.end - location.start + 1) / match.sequenceLength  // Fixed calculation
                double r = Math.abs(location.hmmEnd - location.hmmStart + 1) / (location.envelopeEnd - location.envelopeStart + 1)  // Fixed calculation
                double ld = Math.abs(match.sequenceLength - datEntries[modelAccession].meanL)

                if (datChildren.containsKey(modelAccession)) {
                    if (r > LENGTH_RATIO_THRESHOLD && location.score >= datEntries[modelAccession].minS) {
                        String parent = datChildren[modelAccession]
                        if (store[proteinAccession]?.containsKey(parent)) {
                            processedMatches.computeIfAbsent(proteinAccession, { [:] })[parent] = store[proteinAccession][parent]
                        }
                        else {
                            promote["${proteinAccession}-${parent}"] = true
                        }
                        UpdateMatch(processedMatches, proteinAccession, modelAccession, match, location)
                    }
                    else if (
                        r > LENGTH_RATIO_THRESHOLD &&
                        ovl >= OVERLAP_THRESHOLD &&
                        location.score >= datEntries[modelAccession].minS &&
                        (ld < LENGTH_DEVIATION_THRESHOLD * datEntries[modelAccession].stdL || ld < MINIMUM_LENGTH_DEVIATION)
                    ) {
                        UpdateMatch(processedMatches, proteinAccession, modelAccession, match, location)
                    }
                    else if (promote["${proteinAccession}-${modelAccession}"]) {
                        UpdateMatch(processedMatches, proteinAccession, modelAccession, match, location)
                    }
                    else {
                        store.computeIfAbsent(proteinAccession, { [:] })
                        if (!store[proteinAccession].containsKey(modelAccession)) {
                            store[proteinAccession][modelAccession] = new Match(
                                match.modelAccession,
                                match.evalue,
                                match.score,
                                match.bias
                            )
                        }
                        store[proteinAccession][modelAccession].addLocation(location)
                    }
                }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pirsf.json")
    def json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

def UpdateMatch(Map<String, Map<String, Match>> processedMatches, String proteinAccession, String modelAccession, Match match, Location location) {
    processedMatches.computeIfAbsent(proteinAccession, { [:] })
    if (!processedMatches[proteinAccession].containsKey(modelAccession)) {
        processedMatches[proteinAccession][modelAccession] = new Match(
            match.modelAccession,
            match.evalue,
            match.score,
            match.bias
        )
    }
    processedMatches[proteinAccession][modelAccession].addLocation(location)
}

def processMatchData()
