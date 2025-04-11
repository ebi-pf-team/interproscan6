import groovy.json.JsonOutput

process SEARCH_PIRSF {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.dtbl")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        --domtblout hmmsearch.dtbl \
        ${hmmdb} ${fasta}
    """
}

process PARSE_PIRSF {
    label 'run_locally'

    input:
    tuple val(meta), val(hmmsearch_dtbl)
    val pirsf_dat_file

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def datData = PirsfDatEntry.parsePirsfDatFile(pirsf_dat_file.toString())  // [datEntries, datChildren]
    SignatureLibraryRelease library = new SignatureLibraryRelease("PIRSF", null)

    /* Retrieve the matches from the dtbl file because we need the tlen
    for the hmmLength, and PIRSF does not assign the same hmmer start/ends
    to the same Match.Location fields as other member database. For example,
    PIRSF uses the envelope start and end for the location start/end. */
    Map<String, Map<String, Match>> rawMatches = new LinkedHashMap<>()
    File hmmerDtblFile = new File(hmmsearch_dtbl.toString())
    hmmerDtblFile.withReader { reader ->
        String line
        while ((line = reader.readLine()) != null) {
            if (!line.startsWith("#")) {
                def lineData = line.trim().split("\\s+")
                String proteinAccession = lineData[0]
                String modelAccession = lineData[4]

                rawMatches.computeIfAbsent(proteinAccession, {new LinkedHashMap<>()})
                if (!rawMatches[proteinAccession].containsKey(modelAccession)) {
                    Match match = new Match(
                        modelAccession,
                        lineData[6].toDouble(),  // seqEvalue
                        lineData[7].toDouble(),  // seqScore
                        lineData[8].toDouble()   // seqBias
                    )
                    match.sequenceLength = lineData[5].toInteger()
                    rawMatches[proteinAccession].put(modelAccession, match)
                }

                Location location = new Location(
                    lineData[19].toInteger(),   // Start: use the envelope
                    lineData[20].toInteger(),   // End: use the envelope
                    lineData[17].toInteger(),   // hmmStart: use ali from
                    lineData[18].toInteger(),   // hmmEnd: use ali to
                    lineData[5].toInteger(),    // hmmLength
                    null,                       // hmmBounds -- added later when merging matches together
                    lineData[19].toInteger(),   // envelopeStart
                    lineData[20].toInteger(),   // envelopeEnd
                    lineData[12].toDouble(),    // domain e-value
                    lineData[13].toFloat(),     // domain score
                    lineData[14].toFloat()      // domain bias
                )
                rawMatches[proteinAccession][modelAccession].addLocation(location)
            }
        }
    }

    /* Filter the matches using the data in the dat file */
    def datEntries = datData[0]
    def datChildren = datData[1]

    double LENGTH_RATIO_THRESHOLD = 0.67
    double OVERLAP_THRESHOLD = 0.8
    double LENGTH_DEVIATION_THRESHOLD = 3.5
    int MINIMUM_LENGTH_DEVIATION = 50

    Map<String, Map<String, Match>> processedMatches = new LinkedHashMap<>()
    def promote = [:]
    def store = [:]

    rawMatches.each { proteinAccession, modelMatches ->
        modelMatches.each { modelAccession, rawMatch ->
            // Merge match locations together
            match = processMatchLocations(rawMatch, library)

            // Calculate ratios
            // Overall length
            double ovl = Math.abs(match.locations[0].end - match.locations[0].start + 1) / match.sequenceLength
            // Ratio over coverage of sequence and profile hmm
            double r = Math.abs(match.locations[0].hmmEnd - match.locations[0].hmmStart + 1) / (match.locations[0].envelopeEnd - match.locations[0].envelopeStart + 1)
            // length deviation
            double ld = Math.abs(match.sequenceLength - datEntries[modelAccession].meanL)

            if (datChildren.containsKey(modelAccession)) {
                // process a subfamily match
                if (r > LENGTH_RATIO_THRESHOLD && match.locations[0].score >= datEntries[modelAccession].minS) {
                    String parent = datChildren[modelAccession]
                    if (store[proteinAccession]?.containsKey(parent)) {
                        processedMatches.computeIfAbsent(proteinAccession, { [:] })[parent] = store[proteinAccession][parent]
                    }
                    else {
                        promote["${proteinAccession}-${parent}"] = true
                    }
                    UpdateMatch(processedMatches, proteinAccession, match)
                }
            }

            else if (
                r > LENGTH_RATIO_THRESHOLD &&
                ovl >= OVERLAP_THRESHOLD &&
                match.locations[0].score >= datEntries[modelAccession].minS &&
                (ld < LENGTH_DEVIATION_THRESHOLD * datEntries[modelAccession].stdL || ld < MINIMUM_LENGTH_DEVIATION)
            ) {
                UpdateMatch(processedMatches, proteinAccession, match)
            }

            else if (promote["${proteinAccession}-${modelAccession}"]) {
                UpdateMatch(processedMatches, proteinAccession, match)
            }

            else {
                // store for potential use
                store.computeIfAbsent(proteinAccession, { [:] })
                store[proteinAccession].computeIfAbsent(modelAccession, {match})
            }
        }
    }

    /* filter to use only the best matches, by score */
    def bestMatches = [:]
    processedMatches.each { proteinAccession, proteinMatches ->
        def matchesSorted = proteinMatches.keySet().sort { key ->
            proteinMatches[key]?.locations?.get(0)?.score ?: 0
        }.reverse()

        matchesSorted.each { modelAccession ->
            // ignore subfamilies
            if (modelAccession ==~ /^PIRSF5/) { return }

            bestMatches.computeIfAbsent(proteinAccession, { [:] })
            if (!bestMatches[proteinAccession].containsKey(modelAccession)) {
                bestMatches[proteinAccession][modelAccession] = proteinMatches[modelAccession]
            }

            // see if this model has subfamilies that also matched the protein
            if (datEntries[modelAccession].children != null) {
                datEntries[modelAccession].children.each { subFamAccession ->
                    if (proteinMatches.containsKey(subFamAccession)) {
                        if (!bestMatches[proteinAccession]?.containsKey(subFamAccession)) {
                            bestMatches[proteinAccession] = bestMatches[proteinAccession] ?: [:]
                            bestMatches[proteinAccession][subFamAccession] = proteinMatches[subFamAccession]
                        }
                    }
                }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pirsf.json")
    def json = JsonOutput.toJson(bestMatches)
    new File(outputFilePath.toString()).write(json)
}

def processMatchLocations(Match match, SignatureLibraryRelease library) {
    /* Combine multiple overlapping or related matches into a single consolidated match region,
    but only when both the sequence and HMM model agree on the extensions.
    Sequence:  1....5....10...15...20...25...30...35...40
    Match 1:              |-------|            # seq: 15-25, hmm: 10-20
    Match 2:           |-----|                 # seq: 12-18, hmm: 8-15
    Match 3:                     |----------|  # seq: 22-35, hmm: 18-30
    Combined:          |--------------------|  # seq: 12-35, hmm: 8-30
    */
    // Initialise with the first location
    int seqStart = match.locations[0].start
    int seqEnd = match.locations[0].end
    int hmmStart = match.locations[0].hmmStart
    int hmmEnd = match.locations[0].hmmEnd
    // check subsequent matches
    match.locations.each { location ->
        if (location.start < seqStart && location.hmmStart < hmmStart) {
            seqStart = location.start
            hmmStart = location.hmmStart
        }
        if (location.end > seqEnd && location.hmmEnd < hmmEnd) {
            seqEnd = location.end
            hmmEnd = location.hmmEnd
        }
    }

    Match processedMatch = new Match(
        match.modelAccession,
        match.evalue,
        match.score,
        match.bias,
        new Signature(match.modelAccession, library)
    )
    processedMatch.sequenceLength = match.sequenceLength
    String hmmBoundStart = hmmStart == 1 ? "[" : "."
    String hmmBoundEnd = hmmEnd == match.locations[0].hmmLength ? "]" : "."
    processedMatch.addLocation(
        new Location(
            seqStart,
            seqEnd,   // start + end == envelopeStart + envelopeEnd
            hmmStart,
            hmmEnd,
            match.locations[0].hmmLength,
            "${hmmBoundStart}${hmmBoundEnd}",
            seqStart,
            seqEnd,   // envelopeStart + envelopeEnd
            match.locations[0].evalue,
            match.locations[0].score,
            match.locations[0].bias
        )
    )

    return processedMatch
}

def UpdateMatch(Map<String, Map<String, Match>> processedMatches, String proteinAccession, Match match) {
    processedMatches.computeIfAbsent(proteinAccession, { [:] })
    if (!processedMatches[proteinAccession].containsKey(match.modelAccession)) {
        processedMatches[proteinAccession][match.modelAccession] = match
    }
}
