import groovy.json.JsonOutput

process RUN_HMMER {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmdb
    val options    // e.g. "-Z 65245 -E 0.001"

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        ${options} \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}
