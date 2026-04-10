import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.CATH
import uk.ac.ebi.interpro.HMMER3
import uk.ac.ebi.interpro.Match

process SEARCH_GENE3D {
    label     'mem_min', 'time_short', 'dynamic'
    container 'interpro/cath:0.16.10'

    input:
    tuple val(meta), val(meta2), path(fasta)
    tuple path(hmmfile), path(dom2fam), path(disc_pickle)

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out"), path("cath.tsv")

    script:
    """
    hmmsearch \
        -Z 65245 -E 0.001 \
        --cpu ${task.cpus} \
        ${hmmfile} ${fasta} > hmmsearch.out

    cath-resolve-hits \
        --input-format=hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore=25 \
        --output-hmmer-aln hmmsearch.out > resolved.out

    assign_cath_superfamilies.py \
        ${dom2fam} \
        ${disc_pickle} \
        resolved.out \
        cath.tsv
    """
}

process PARSE_CATHGENE3D {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(hmmseach_out), val(cath_tsv)

    output:
    tuple val(meta), val(meta2), path("cathgene3d.json")

    exec:
    def member_db = "CATH-Gene3D"
    def hmmer_matches = HMMER3.parseOutput(hmmseach_out, member_db)
    def cath_domains = CATH.parseAssignedFile(cath_tsv)
    def matches = CATH.mergeWithHmmerMatches(cath_domains, hmmer_matches, member_db)
    def filepath = task.workDir.resolve("cathgene3d.json")
    filepath.text = JsonOutput.toJson(matches)
}

process PREPARE_FUNFAM {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(cathgene3d_json)

    output:
    tuple val(meta), val(meta2), path("cath.txt")
    
    exec:
    def seen = [] as Set
    task.workDir
        .resolve("cath.txt")
        .withWriter { writer ->
            new JsonSlurper().parse(cathgene3d_json)
                .values()
                .each { jsonMatches ->
                    jsonMatches
                        .values()
                        .each { jsonMatch ->
                            def match = Match.fromMap(jsonMatch)
                            def accession = match?.signature?.accession
                            assert accession != null
                            assert accession.startsWith("G3DSA:")
                            def cath_id = accession.substring(6)
                            if (seen.add(cath_id)) {
                                writer.writeLine(cath_id)
                            }
                        }
                }
        }
}

process SEARCH_FUNFAM {
    label     'mem_min', 'time_short', 'dynamic'
    container 'interpro/cath:0.16.10'

    input:
    tuple val(meta), val(meta2), path(fasta), path(caths_txt)
    path root_dir

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out"), path("resolved.out")

    script:
    """
    touch hmmsearch.out
    while IFS= read -r id || [[ -n "\$id" ]]; do
        [[ -z "\$id" ]] && continue
        
        hmm="${root_dir}/\$(tr '.' '/' <<< "\$id").hmm"
        
        if [[ -f "\$hmm" ]]; then
            hmmsearch -Z 65245 --cut_tc --cpu ${task.cpus} "\$hmm" ${fasta} >> hmmsearch.out
        fi
    done < ${caths_txt}

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
