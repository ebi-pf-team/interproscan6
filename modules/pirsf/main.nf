import groovy.json.JsonOutput

process SEARCH_PIRSF {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dir
    val hmm

    output:
    tuple val(meta), path("hmmsearch.dtbl")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        --domtblout hmmsearch.dtbl \
        ${dir}/${hmm} ${fasta}
    """
}

process PARSE_PIRSF {
    executor 'local'

    input:
    tuple val(meta), val(hmmsearch_dtbl)
    val dirpath
    val datfile

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def datPath = "${dirpath.toString()}/${datfile}"
    def datData = PirsfDatEntry.parsePirsfDatFile(datPath)  // [datEntries, datChildren]
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
                    match.sequenceLength = lineData[2].toInteger()
                    rawMatches[proteinAccession].put(modelAccession, match)
                }
                Location location = new Location(
                    lineData[17].toInteger(),   // Start
                    lineData[18].toInteger(),   // End
                    lineData[15].toInteger(),   // hmmStart
                    lineData[16].toInteger(),   // hmmEnd
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

            /* Combine multiple overlapping or related matches into a single consolidated match region,
            but only when both the sequence and HMM model agree on the extensions.
            Sequence:  1....5....10...15...20...25...30...35...40
            Match 1:              |-------|            # seq: 15-25, hmm: 10-20
            Match 2:           |-----|                 # seq: 12-18, hmm: 8-15
            Match 3:                     |----------|  # seq: 22-35, hmm: 18-30
            Combined:          |--------------------|  # seq: 12-35, hmm: 8-30
            */
            int seqStart = rawMatch.locations[0].start
            int seqEnd = rawMatch.locations[0].end
            int hmmStart = rawMatch.locations[0].hmmStart
            int hmmEnd = rawMatch.locations[0].hmmEnd

            rawMatch.locations.each { location ->
                if (location.start < seqStart && location.hmmStart < hmmStart) {
                    seqStart = location.start
                    hmmStart = location.hmmStart
                }
                if (location.end > seqEnd && location.hmmEnd > hmmEnd) {
                    seqEnd = location.end
                    hmmEnd = location.hmmEnd
                }
            }

            // Calculate ratios
            // Overall length
            double ovl = Math.abs(seqEnd - seqStart + 1) / rawMatch.sequenceLength
            // Ratio over coverage of sequence and profile hmm
            double r = Math.abs(hmmEnd - hmmStart + 1) / (seqEnd - seqStart + 1)
            // length deviation
            double ld = Math.abs(rawMatch.sequenceLength - datEntries[modelAccession].meanL)

            match = processMatchLocations(rawMatch, library)

            if (datChildren.containsKey(modelAccession)) {
                // process a subfamily match
                if (r > LENGTH_RATIO_THRESHOLD && match.score >= datEntries[modelAccession].minS) {
                    String parent = datChildren[modelAccession]
                    if (store[proteinAccession]?.containsKey(parent)) {
                        processedMatches.computeIfAbsent(proteinAccession, { [:] })[parent] = store[proteinAccession][parent]
                    } else {
                        promote["${proteinAccession}-${parent}"] = true
                    }
                    UpdateMatch(processedMatches, proteinAccession, match)
                }
            } else if (
                r > LENGTH_RATIO_THRESHOLD &&
                    ovl >= OVERLAP_THRESHOLD &&
                    match.score >= datEntries[modelAccession].minS &&
                    (ld < LENGTH_DEVIATION_THRESHOLD * datEntries[modelAccession].stdL || ld < MINIMUM_LENGTH_DEVIATION)
            ) {
                UpdateMatch(processedMatches, proteinAccession, match)
            } else if (promote["${proteinAccession}-${modelAccession}"]) {
                UpdateMatch(processedMatches, proteinAccession, match)
            } else {
                // store for potential use
                store.computeIfAbsent(proteinAccession, { [:] })
                store[proteinAccession].computeIfAbsent(modelAccession, {match})
            }
        }
    }

    /* filter to use only the best matches, by evalue ASC, score DESC */
    def bestMatches = [:]
    processedMatches.each { proteinAccession, proteinMatches ->
        def matchesSorted = proteinMatches.keySet().sort { key ->
            proteinMatches[key].score ?: 0
        }.reverse()

        // Find the first non-subfamily match (the best one)
        def bestMatch = matchesSorted.find { !(it ==~ /^PIRSF5/) }
        bestMatches.computeIfAbsent(proteinAccession, { [:] })
        bestMatches[proteinAccession][bestMatch] = proteinMatches[bestMatch]

        // see if this model has subfamilies that also matched the protein
        if (datEntries[bestMatch]?.children != null) {
            datEntries[bestMatch].children.each { subFamAccession ->
                if (proteinMatches.containsKey(subFamAccession)) {
                    bestMatches[proteinAccession][subFamAccession] = proteinMatches[subFamAccession]
                }
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pirsf.json")
    def json = JsonOutput.toJson(bestMatches)
    new File(outputFilePath.toString()).write(json)
}

def processMatchLocations(Match match, SignatureLibraryRelease library) {
    Match processedMatch = new Match(
        match.modelAccession,
        match.evalue,
        match.score,
        match.bias,
        new Signature(match.modelAccession, library)
    )
    processedMatch.sequenceLength = match.sequenceLength
    match.locations.each { location ->
        String hmmBoundStart = location.hmmStart == 1 ? "[" : "."
        String hmmBoundEnd = location.hmmEnd == location.hmmLength ? "]" : "."
        processedMatch.addLocation(
            new Location(
                location.envelopeStart,  // start is the envelopeStart value on output
                location.envelopeEnd,   // end is the envelopeEnd value on output
                location.start,  // hmmStart is the start value on output
                location.end,  // hmmEnd is the end value on output
                location.hmmLength,
                "${hmmBoundStart}${hmmBoundEnd}",
                location.envelopeStart,
                location.envelopeEnd,
                location.evalue,
                location.score,
                location.bias
            )
        )
    }
    return processedMatch
}

def UpdateMatch(Map<String, Map<String, Match>> processedMatches, String proteinAccession, Match match) {
    processedMatches.computeIfAbsent(proteinAccession, { [:] })
    if (!processedMatches[proteinAccession].containsKey(match.modelAccession)) {
        processedMatches[proteinAccession][match.modelAccession] = match
    }
}
