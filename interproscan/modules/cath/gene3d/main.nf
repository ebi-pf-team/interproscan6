import groovy.json.JsonOutput

process SEARCH_GENE3D {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 65245 -E 0.001 \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process RESOLVE_GENE3D {
    label 'analysis_parser'

    input:
    tuple val(meta), path(hmmseach_out)

    output:
    tuple val(meta), path("resolved.out")

    // cath_resolve_hits is a third party tool used to minimise suprious hits
    script:
    """
    /opt/cath-tools/cath-resolve-hits \
        ${hmmseach_out} \
        --input-for hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore 25 \
        --output-hmmer-aln > resolved.out
    """
}

process ASSIGN_CATH_SUPFAM {
    label 'analysis_parser'

    input:
    tuple val(meta), path(cath_resolve_out)
    path dom2fam_file
    path disc_pickle_file

    output:
    tuple val(meta), path("cath.tsv")

    script:
    """
    python ${projectDir}/bin/cath/assign_cath_superfamilies.py \
        ${dom2fam_file} \
        ${disc_pickle_file} \
        ${cath_resolve_out} \
        cath.tsv
    """
}

process PARSE_CATHGENE3D {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmseach_out), val(cath_tsv)

    output:
    tuple val(meta), path("cathgene3d.json")

    exec:
    def gene3dMatches = HMMER3.parseOutput(hmmseach_out.toString())

    def gene3dDomains = [:]
    gene3dMatches.each { sequenceId, matches ->
        matches.values().each { m1 ->
            m1.locations.each { loc ->
                Match m2 = new Match(
                    m1.modelAccession, 
                    m1.evalue,
                    m1.score, 
                    m1.bias)
                m2.addLocation(loc)
                String key = "${m2.modelAccession}-${loc.envelopeStart}-${loc.envelopeEnd}"

                if (gene3dDomains.containsKey(sequenceId)) {
                    assert !gene3dDomains.containsKey(key)
                    gene3dDomains[sequenceId][key] = m2
                } else {
                    gene3dDomains[sequenceId] = [key: m2]
                }
            }
        }
    }

    def cathDomains = [:]
    cath_tsv.eachLine { line ->
        if (line[0] != "#") {
            def fields = line.split("\t")
            assert fields.size() == 10
            String sequenceId = fields[2]
            def dom = new CathDomain(
                fields[0],                                                  // Domain ID
                fields[3],                                                  // Match ID
                fields[1],                                                  // CATH ID
                Double.parseDouble(fields[4]),                              // Score
                Double.parseDouble(fields[9]),                              // i-evalue
                fields[5].split(",").collect { new SimpleLocation(it) },    // Boundaries
                fields[6].split(",").collect { new SimpleLocation(it) },    // Resolved
            )

            if (cathDomains.containsKey(sequenceId)) {
                cathDomains[sequenceId].add(dom)
            } else {
                cathDomains[sequenceId] = [dom]
            }
        }
    }

    def matches = [:]
    cathDomains.each { sequenceId, domains ->
        def hmmerDomains = gene3dDomains.get(sequenceId)
        assert hmmerDomains != null

        domains.each { cathDomain ->
            def key = cathDomain.getKey()
            def hmmerDomain = hmmerDomains.get(key)
            assert hmmerDomain != null

            def boundaries = cathDomain.resolvedBoundaries
            def fragments = []
            if (boundaries.size() > 1) {
                cathDomain.resolvedBoundaries.eachWithIndex { x, i ->
                    LocationFragment fragment
                    if (i == 0) {
                        fragment = new LocationFragment(x.start, x.end, "C_TERMINAL_DISC")
                    } else if (i + 1 < boundaries.size()) {
                        fragment = new LocationFragment(x.start, x.end, "NC_TERMINAL_DISC")
                    } else {
                        fragment = new LocationFragment(x.start, x.end, "N_TERMINAL_DISC")
                    }
                    fragments.add(fragment)
                }
            } else {
                LocationFragment fragment = new LocationFragment(boundaries[0].start, boundaries[0].end, "CONTINUOUS")
                fragments.add(fragment)
            }

            Location location = new Location(
                cathDomain.getStart(),
                cathDomain.getEnd(),
                hmmerDomain.locations[0].hmmStart,
                hmmerDomain.locations[0].hmmEnd,
                hmmerDomain.locations[0].hmmLength,
                hmmerDomain.locations[0].hmmBounds,
                hmmerDomain.locations[0].envelopeStart,
                hmmerDomain.locations[0].envelopeEnd,
                cathDomain.evalue,
                cathDomain.score,
                null,
                fragments
            )

            // domain.addLocation(location)

            def sequenceDomains
            if (matches.containsKey(sequenceId)) {
                sequenceDomains = matches[sequenceId]
            } else {
                sequenceDomains = (matches[sequenceId] = [:]) 
            }

            def domId = hmmerDomain.modelAccession
            if (sequenceDomains.containsKey(domId)) {
                sequenceDomains[domId].addLocation(location)
            } else {
                Match domain = new Match(
                    cathDomain.domainId, 
                    hmmerDomain.evalue,
                    hmmerDomain.score, 
                    hmmerDomain.bias
                )
                domain.signature = new Signature(cathDomain.supfamId)
                domain.addLocation(location)
                sequenceDomains[domId] = domain
            }
        }
    }

    def outputFilePath = task.workDir.resolve("cathgene3d.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}