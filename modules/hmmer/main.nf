process RUN_HMMER {
    label     'mem_min', 'time_medium', 'dynamic'
    container 'interpro/hmmer:3.3'

    input:
    tuple val(meta), path(fasta)
    path hmm
    val options    // e.g. "-Z 65245 -E 0.001"

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    hmmsearch \
        ${options} \
        --cpu ${task.cpus} \
        ${hmm} ${fasta} > hmmsearch.out
    """
}
