import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.CATH
import uk.ac.ebi.interpro.HMMER3
import uk.ac.ebi.interpro.Match

process PREPARE_FUNFAM {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(cathgene3d_json)

    output:
    tuple val(meta), val(meta2), val(funfams)
    
    exec:
    funfams = new JsonSlurper().parse(cathgene3d_json)
        .values()
        .collect{ jsonMatches ->
            jsonMatches
                .values()
                .collect { jsonMatch ->
                    def match = Match.fromMap(jsonMatch)
                    def accession = match?.signature?.accession                       
                    assert accession != null
                    assert accession.startsWith("G3DSA:")
                    def cathId = accession.substring(6)
                    return cathId
                }
        }
        .flatten()
        .unique()
}

process SEARCH_FUNFAM {
    label 'mem_min', 'time_short', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(fasta), val(funfams)
    path root_dir

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out"), path("resolved.out")

    script:
    def commands = "touch hmmsearch.out\n"
    funfams.each { cathId -> 
        def leaf = cathId
            .split("\\.")  // codenarc-disable-line JoinMismatchRule, JoinDuplicateRule
            .join(File.separator) + ".hmm"
        def hmm = root_dir.resolve(leaf)
        if (hmm.exists()) {
            commands += "hmmsearch"
            commands += " -Z 65245 --cut_tc"
            commands += " --cpu ${task.cpus}"
            commands += " ${hmm} ${fasta} >> hmmsearch.out\n"
        }
    }

    """
    ${commands}

    cath-resolve-hits \
        --input-format=hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore=25 \
        --output-hmmer-aln hmmsearch.out > resolved.out
    """
}


process PARSE_FUNFAM {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_out), val(resolved_tsv)

    output:
    tuple val(meta), val(meta2), path("cathfunfam.json")

    exec:
    def member_db = "CATH-FunFam"
    def hmmer_matches = HMMER3.parseOutput(hmmseach_out, member_db)
    def funfam_domains = CATH.parseResolvedFile(resolved_tsv)
    def matches = CATH.mergeWithHmmerMatches(funfam_domains, hmmer_matches, member_db)
    def filepath = task.workDir.resolve("cathfunfam.json")
    filepath.text = JsonOutput.toJson(matches)
}