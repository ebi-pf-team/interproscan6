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
    SignatureLibraryRelease library = new SignatureLibraryRelease("PIRSF", null)

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
        print(sequences[proteinAccession])
        sequenceLength = sequences[proteinAccession].length()
        modelMatches.each { modelAccession, rawMatch ->
            int seqStart = Integer.MAX_VALUE
            int seqEnd = Integer.MIN_VALUE
            int hmmStart = Integer.MAX_VALUE
            int hmmEnd = Integer.MIN_VALUE
            rawMatch.locations.each { location ->
                if (location.included) {
                    if (location.start < seqStart && location.hmmStart < hmmStart) {
                        seqStart = location.start
                        hmmStart = location.hmmStart
                    }
                    if (location.end > seqEnd && location.hmmEnd > hmmEnd) {
                        seqEnd = location.end
                        hmmEnd = location.hmmEnd
                    }
                }
            }
            if (seqStart == Integer.MAX_VALUE && seqEnd == Integer.MIN_VALUE) {
                // No valid locations found
                return
            }

            // Calculate ratios
            // Overall length
            double ovl = Math.abs(seqEnd - seqStart + 1) / sequenceLength
            // Ratio over coverage of sequence and profile hmm
            double r = Math.abs(hmmEnd - hmmStart + 1) / (seqEnd - seqStart + 1)
            // length deviation
            double ld = Math.abs(sequenceLength - datEntries[modelAccession].meanL)
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
                location.start,
                location.end,
                location.hmmStart,
                location.hmmEnd,
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
