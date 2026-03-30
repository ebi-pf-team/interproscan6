import groovy.json.JsonOutput
import uk.ac.ebi.interpro.CATH
import uk.ac.ebi.interpro.HMMER3
import uk.ac.ebi.interpro.Match

process SEARCH_GENE3D {
    label 'mem_min', 'time_short', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), val(meta2), path(fasta)
    path dirpath
    val hmmfile
    val dom2fam
    val disc_pickle

    output:
    tuple val(meta), val(meta2), path("hmmsearch.out"), path("cath.tsv")

    script:
    """
    hmmsearch \
        -Z 65245 -E 0.001 \
        --cpu ${task.cpus} \
        ${dirpath}/${hmmfile} ${fasta} > hmmsearch.out

    cath-resolve-hits \
        --input-format=hmmsearch_out \
        --min-dc-hmm-coverage=80 \
        --worst-permissible-bitscore=25 \
        --output-hmmer-aln hmmsearch.out > resolved.out

    python ${projectDir}/bin/cath/assign_cath_superfamilies.py \
        ${dirpath}/${dom2fam} \
        ${dirpath}/${disc_pickle} \
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