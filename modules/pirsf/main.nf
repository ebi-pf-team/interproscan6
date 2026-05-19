process SEARCH_PIRSF {
    label     'mem_low', 'time_short', 'dynamic'
    container 'interpro/hmmer:3.3'

    input:
    tuple val(meta), path(fasta)
    path hmm

    output:
    tuple val(meta), path(fasta), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        -E 0.01 --acc \
        --cpu ${task.cpus} \
        ${hmm} ${fasta} > hmmsearch.out
    """
}

process PARSE_PIRSF {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(fasta), val(hmmsearch_out)
    val datfile

    output:
    tuple val(meta), path("pirsf.json")

    exec:
    def (models, subfamilies) = parseDatFile(datfile)
    def hmmerMatches = uk.ac.ebi.interpro.HMMER3.parseOutput(hmmsearch_out, "PIRSF")
    def sequences = uk.ac.ebi.interpro.FastaFile.parse(fasta)

    def LENGTH_RATIO_THRESHOLD = 0.67
    def OVERLAP_THRESHOLD = 0.8
    def LENGTH_DEVIATION_THRESHOLD = 3.5
    def MINIMUM_LENGTH_DEVIATION = 50

    def results = [:] // proteinAccession -> modelAccession -> Match

    hmmerMatches.each { proteinAccession, modelMatches ->
        def sequenceLength = sequences[proteinAccession].length()

        def familyMatches = []
        def subfamilyMatches = []

        modelMatches.each { modelAccession, rawMatch ->
            def seqStart = Integer.MAX_VALUE
            def seqEnd = Integer.MIN_VALUE
            def hmmStart = Integer.MAX_VALUE
            def hmmEnd = Integer.MIN_VALUE
            def envStart = 0
            def envEnd = 0
            def locationScore = 0.0
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
            def (meanL, stdL, minS, _meanS, _stdS) = models[match.modelAccession]

            def location = match.locations[0]
            
            // Overall length
            def ovl = (location.end - location.start + 1) / sequenceLength

            // Length deviation
            def ld = Math.abs(sequenceLength - meanL)

            // Ratio over coverage of sequence and profile hmm
            def r = (location.hmmEnd - location.hmmStart + 1) / (location.end - location.start + 1)

            return r > LENGTH_RATIO_THRESHOLD
                && ovl >= OVERLAP_THRESHOLD
                && match.score >= minS
                && (ld < LENGTH_DEVIATION_THRESHOLD * stdL || ld < MINIMUM_LENGTH_DEVIATION)
        }

        // Select best family match
        def familyMatch = filteredFamilyMatches ? filteredFamilyMatches.max { m -> m.score } : null

        // Filter subfamily matches and select the best one
        def filteredSubfamilyMatches = subfamilyMatches.findAll { match ->
            // Only consider subfamily matches that have a parent family matching the sequence
            def parentFamily = subfamilies[match.modelAccession]

            if (familyMatch && familyMatch.modelAccession != parentFamily) {
                // If we already have a best family match, we want to subfamily to belong the this family
                return false
            }

            def parentMatch = familyMatches.find { m -> m.modelAccession == parentFamily }
            if (!parentMatch) {
                // We don't want a subfamily match if the parent family doesn't hit the sequence
                return false
            }

            def (_meanL, _stdL, minS, _meanS, _stdS) = models[match.modelAccession]
            
            // Ratio over coverage of sequence and profile hmm
            def location = match.locations[0]
            def r = (location.hmmEnd - location.hmmStart + 1) / (location.end - location.start + 1)
            return r > LENGTH_RATIO_THRESHOLD && match.score >= minS
        }

        def subfamilyMatch = filteredSubfamilyMatches ? filteredSubfamilyMatches.max { m -> m.score } : null
        
        if (subfamilyMatch && !familyMatch) {
            // Promote parent family (even if it didn't pass the cutoffs)
            def parentFamily = subfamilies[subfamilyMatch.modelAccession]
            familyMatch = familyMatches.find { m -> m.modelAccession == parentFamily }
        }

        if (familyMatch) {
            results[proteinAccession] = [
                (familyMatch.modelAccession): familyMatch
            ]

            if (subfamilyMatch) {
                results[proteinAccession][subfamilyMatch.modelAccession] = subfamilyMatch
            }
        }
    }

    def filepath = task.workDir.resolve("pirsf.json")
    filepath.text  = groovy.json.JsonOutput.toJson(results)
}

def createMatch(
    match,
    seqStart,
    seqEnd,
    hmmStart,
    hmmEnd,
    envStart,
    envEnd,
    locationScore) {
    def processedMatch = new uk.ac.ebi.interpro.Match(
        match.modelAccession,
        match.evalue,
        locationScore,
        match.bias,
        match.signature
    )
    def hmmBoundStart = hmmStart == 1 ? "[" : "."
    def hmmBoundEnd = hmmEnd == match.locations[0].hmmLength ? "]" : "."
    processedMatch.addLocation(
        new uk.ac.ebi.interpro.Location(
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

def parseDatFile(datFile) {
    def models = [:]    // PIRSF -> list of 5 doubles (meanL, stdL, minS, meanS, stdS)
    def subfamilies = [:]   // child PIRSF -> parent PIRSF

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