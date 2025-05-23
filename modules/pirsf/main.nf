import groovy.json.JsonOutput

process SEARCH_PIRSF {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dir
    val hmm

    output:
    tuple val(meta), path(fasta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        ${dir}/${hmm} ${fasta} > hmmsearch.out
    """
}

process PARSE_PIRSF {
    executor 'local'

    input:
    tuple val(meta), val(fasta), val(hmmsearch_out)
    val dirpath
    val datfile

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def datPath = "${dirpath.toString()}/${datfile}"
    def datData = PirsfDatEntry.parsePirsfDatFile(datPath)  // [datEntries, datChildren]
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "PIRSF")
    Map<String, String> sequences = FastaFile.parse(fasta.toString())

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

    hmmerMatches.each { proteinAccession, modelMatches ->
        int sequenceLength = sequences[proteinAccession].length()
        modelMatches.each { modelAccession, rawMatch ->
            int seqStart = Integer.MAX_VALUE
            int seqEnd = Integer.MIN_VALUE
            int hmmStart = Integer.MAX_VALUE
            int hmmEnd = Integer.MIN_VALUE
            int envStart = 0
            int envEnd = 0
            rawMatch.locations.each { location ->
                if (location.included) {
                    if (location.start < seqStart && location.hmmStart < hmmStart) {
                        seqStart = location.start
                        hmmStart = location.hmmStart
                        envStart = location.envelopeStart
                    }
                    if (location.end > seqEnd && location.hmmEnd > hmmEnd) {
                        seqEnd = location.end
                        hmmEnd = location.hmmEnd
                        envEnd = location.envelopeEnd
                    }
                }
            }
            if (seqStart == Integer.MAX_VALUE && seqEnd == Integer.MIN_VALUE) {
                // No valid locations found
                return
            }

            // Calculate ratios
            // Overall length
            double ovl = (seqEnd - seqStart + 1) / sequenceLength
            // Ratio over coverage of sequence and profile hmm
            double r = Math.abs(hmmEnd - hmmStart + 1) / (seqEnd - seqStart + 1)
            // length deviation
            double ld = Math.abs(sequenceLength - datEntries[modelAccession].meanL)
            match = processMatchLocations(rawMatch, seqStart, seqEnd, hmmStart, hmmEnd)
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

    /* filter to use only the best matches */
    def bestMatches = [:]
    processedMatches.each { proteinAccession, proteinMatches ->
        def matchesSorted = proteinMatches.keySet().sort { key ->
            proteinMatches[key].score ?: 0
        }.reverse()

        // Find the first non-subfamily match (the best one)
        def bestMatch = matchesSorted.find { !it.startsWith("PIRSF5") }
        bestMatches.computeIfAbsent(proteinAccession, { [:] })
        bestMatches[proteinAccession][bestMatch] = proteinMatches[bestMatch]

        // see if this model has subfamilies that also matched the protein
        if (datEntries[bestMatch]?.children) {
            def bestSubfamily = matchesSorted.find { datEntries[bestMatch]?.children.contains(it) }
            if (bestSubfamily) {
                bestMatches[proteinAccession][bestSubfamily] = proteinMatches[bestSubfamily]
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pirsf.json")
    def json = JsonOutput.toJson(bestMatches)
    new File(outputFilePath.toString()).write(json)
}

def processMatchLocations(
        Match match,
        int seqStart,
        int seqEnd,
        int hmmStart,
        int hmmEnd) {
    Match processedMatch = new Match(
        match.modelAccession,
        match.evalue,
        match.score,
        match.bias,
        match.signature
    )
    processedMatch.sequenceLength = match.sequenceLength
    String hmmBoundStart = hmmStart == 1 ? "[" : "."
    String hmmBoundEnd = hmmEnd == match.locations[0].hmmLength ? "]" : "."
    def minEnvelopeStart = match.locations.collect { it.envelopeStart }.min()
    def maxEnvelopeEnd = match.locations.collect { it.envelopeEnd }.max()
    processedMatch.addLocation(
        new Location(
            seqStart,
            seqEnd,
            hmmStart,
            hmmEnd,
            match.locations[0].hmmLength,
            "${hmmBoundStart}${hmmBoundEnd}",
            minEnvelopeStart,
            maxEnvelopeEnd,
            match.evalue,
            match.score,
            match.bias
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
