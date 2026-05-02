process SEARCH_PANTHER {
    label     'mem_low', 'time_short', 'dynamic'
    container 'interpro/panther:1.0'

    input:
    tuple val(meta), val(meta2), path(fasta)
    path hmm

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        -Z 65000000 -E 0.001 \
        --cpu ${task.cpus} \
        ${hmm} ${fasta} > hmmsearch.out
    """
}

process PREPARE_TREEGRAFTER {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_out)

    output:
    tuple val(meta), val(meta2), path("panther.json")
    
    exec:
    def hmmer_matches = uk.ac.ebi.interpro.HMMER3.parseOutput(hmmseach_out, "PANTHER")

    hmmer_matches = hmmer_matches.collectEntries { seqId, matches ->
        // Filter matches to only those with locations that have a score > 100
        def filteredMatches = matches
            .values()
            .collect { m1 ->
                if (!m1.included) {
                    return null
                }

                def locations = m1.locations.findAll { it.included }
                if (!locations) {
                    return null
                }

                // Rename model accession (PTHR23076.orig.30.pir -> PTHR23076)
                def familyId = (m1.modelAccession =~ /^(PTHR\d+)/)[0][1]
                m1.signature.accession = familyId
                def m2 = new uk.ac.ebi.interpro.Match(familyId, m1.evalue, m1.score, m1.bias, m1.signature)
                m2.included = m1.included
                // Only keep the domain with the highest score
                m2.locations = [locations.max { it.score }]
                // Init empty TreeGrafter attribute
                m2.treegrafter = new uk.ac.ebi.interpro.TreeGrafter(null)
                return m2
            }
            .findAll { it != null }

        def bestMatch = filteredMatches.max { it.score }
        return bestMatch ? [(seqId): [(bestMatch.modelAccession): bestMatch]] : [:]
    }

    def jsonFile = task.workDir.resolve("panther.json")
    jsonFile.text = groovy.json.JsonOutput.toJson(hmmer_matches)
}


process RUN_TREEGRAFTER {
    label     'mem_medium', 'time_short', 'dynamic'
    container 'interpro/panther:1.0'
    
    input:
    tuple val(meta), val(meta2), path(json)
    path msf

    output:
    tuple val(meta), val(meta2), path("epang.tsv")

    script:
    """
    treegrafter.py --threads ${task.cpus} ${json} ${msf} > epang.tsv
    """
}

process PARSE_PANTHER {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_json), val(epang_tsv)

    output:
    tuple val(meta), val(meta2), path("panther.json")

    exec:
    def matches = new groovy.json.JsonSlurper().parse(hmmseach_json).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            [(matchId): uk.ac.ebi.interpro.Match.fromMap(jsonMatch)]
        }]
    }

    epang_tsv.eachLine { line ->
        line = line.trim()
        def fields = line.split()
        assert fields.size() == 3
        def seqId = fields[0]
        def familyId = fields[1]
        def nodeId = fields[2]
        
        def match = matches[seqId]?.get(familyId)
        assert match != null
        match.treegrafter = new uk.ac.ebi.interpro.TreeGrafter(nodeId)
    }

    def filepath = task.workDir.resolve("panther.json")
    filepath.text = groovy.json.JsonOutput.toJson(matches)
}