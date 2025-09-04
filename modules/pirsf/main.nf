import groovy.json.JsonOutput

process SEARCH_PIRSF {
    label 'mini', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dir
    val hmm

    output:
    tuple val(meta), path(fasta), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        ${dir}/${hmm} ${fasta} > hmmsearch.out
    """
}

process PARSE_PIRSF {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(fasta), val(hmmsearch_out)
    val dirpath
    val datfile

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def datPath = "${dirpath.toString()}/${datfile}"
    def (models, subfamilies) = parseDatFile(datPath)
    def hmmerMatches = HMMER3.parseOutput(hmmsearch_out.toString(), "PIRSF")
    Map<String, String> sequences = FastaFile.parse(fasta.toString())

    double LENGTH_RATIO_THRESHOLD = 0.67
    double OVERLAP_THRESHOLD = 0.8
    double LENGTH_DEVIATION_THRESHOLD = 3.5
    int MINIMUM_LENGTH_DEVIATION = 50

    def results = [:] // proteinAccession -> modelAccession -> Match

    hmmerMatches.each { proteinAccession, modelMatches ->
        int sequenceLength = sequences[proteinAccession].length()

        def familyMatches = []
        def subfamilyMatches = []

        modelMatches.each { modelAccession, rawMatch ->
            int seqStart = Integer.MAX_VALUE
            int seqEnd = Integer.MIN_VALUE
            int hmmStart = Integer.MAX_VALUE
            int hmmEnd = Integer.MIN_VALUE
            int envStart = 0
            int envEnd = 0
            double locationScore = 0.0
            rawMatch.locations.each { location ->
                if (location.included) {
                    locationScore += location.score
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

            def match = createMatch(rawMatch, seqStart, seqEnd, hmmStart, hmmEnd, envStart, envEnd, locationScore)

            if (subfamilies.containsKey(modelAccession)) {
                subfamilyMatches.add(match)
            } else {
                familyMatches.add(match)
            }
        }

        // Filter family matches
        def filteredFamilyMatches = familyMatches.findAll { match ->
            def (meanL, stdL, minS, meanS, stdS) = models[match.modelAccession]

            Location location = match.locations[0]
            
            // Overall length
            double ovl = (location.end - location.start + 1) / sequenceLength

            // Length deviation
            double ld = Math.abs(sequenceLength - meanL)

            // Ratio over coverage of sequence and profile hmm
            double r = (location.hmmEnd - location.hmmStart + 1) / (location.end - location.start + 1)

            return r > LENGTH_RATIO_THRESHOLD
                && ovl >= OVERLAP_THRESHOLD
                && match.score >= minS
                && (ld < LENGTH_DEVIATION_THRESHOLD * stdL || ld < MINIMUM_LENGTH_DEVIATION)
        }

        // Select best family match
        def familyMatch = filteredFamilyMatches ? filteredFamilyMatches.max { it.score } : null

        // Filter subfamily matches and select the best one
        def filteredSubfamilyMatches = subfamilyMatches.findAll { match ->
            // Only consider subfamily matches that have a parent family matching the sequence
            def parentFamily = subfamilies[match.modelAccession]

            if (familyMatch && familyMatch.modelAccession != parentFamily) {
                // If we already have a best family match, we want to subfamily to belong the this family
                return false
            }

            def parentMatch = familyMatches.find { it.modelAccession == parentFamily }
            if (!parentMatch) {
                // We don't want a subfamily match if the parent family doesn't hit the sequence
                return false
            }

            def (meanL, stdL, minS, meanS, stdS) = models[match.modelAccession]
            
            // Ratio over coverage of sequence and profile hmm
            Location location = match.locations[0]
            double r = (location.hmmEnd - location.hmmStart + 1) / (location.end - location.start + 1)
            return r > LENGTH_RATIO_THRESHOLD && match.score >= minS
        }

        def subfamilyMatch = filteredSubfamilyMatches ? filteredSubfamilyMatches.max { it.score } : null
        
        if (subfamilyMatch && !familyMatch) {
            // Promote parent family (even if it didn't pass the cutoffs)
            def parentFamily = subfamilies[subfamilyMatch.modelAccession]
            familyMatch = familyMatches.find { it.modelAccession == parentFamily }
        }

        if (familyMatch) {
            results[proteinAccession] = [:].with {
                it[familyMatch.modelAccession] = familyMatch
                if (subfamilyMatch) {
                    it[subfamilyMatch.modelAccession] = subfamilyMatch
                }
                it
            }
        }
    }

    def outputFilePath = task.workDir.resolve("pirsf.json")
    def json = JsonOutput.toJson(results)
    new File(outputFilePath.toString()).write(json)
}

def createMatch(
    Match match,
    int seqStart,
    int seqEnd,
    int hmmStart,
    int hmmEnd,
    int envStart,
    int envEnd,
    double locationScore) {
    Match processedMatch = new Match(
        match.modelAccession,
        match.evalue,
        locationScore,
        match.bias,
        match.signature
    )
    String hmmBoundStart = hmmStart == 1 ? "[" : "."
    String hmmBoundEnd = hmmEnd == match.locations[0].hmmLength ? "]" : "."
    processedMatch.addLocation(
        new Location(
            seqStart,
            seqEnd,
            hmmStart,
            hmmEnd,
            match.locations[0].hmmLength,
            "${hmmBoundStart}${hmmBoundEnd}",
            envStart,
            envEnd,
            match.evalue,
            locationScore,
            match.bias
        )
    )
    return processedMatch
}

def parseDatFile(String datPath) {
    def models = [:]    // PIRSF -> list of 5 doubles (meanL, stdL, minS, meanS, stdS)
    def subfamilies = [:]   // child PIRSF -> parent PIRSF

    File datFile = new File(datPath)
    datFile.withReader { reader ->
        def accession = null
        reader.eachLine { line ->
            if (line.startsWith('>')) {
                def parts = line.split(/\s+/)
                accession = parts[0].replace('>', '')

                def match = (line =~ /^>PIRSF\d+\schild:\s(.+)$/)
                if (match) {
                    match[0][1].trim().split(/\s+/).each { child ->
                        subfamilies[child] = accession
                    }
                }

            } else if (line ==~ /\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*/) {
                def values = line.split(/\s+/)*.toDouble()
                models[accession] = values
            }
        }
    }

    return [models, subfamilies]
}
