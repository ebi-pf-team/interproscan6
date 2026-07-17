process RUN_NAIL {
    label     'mem_low'
    label     'time_medium'
    label     'dynamic'
    container 'interpro/nail:0.6.0'

    input:
    tuple val(meta), path(fasta)
    path hmm
    val options    // e.g. "-Z 65245"

    output:
    tuple val(meta), path("nail.tbl")

    script:
    """
    nail search \
        ${options} \
        -t ${task.cpus} \
        --tmp-dir tmp \
        --tbl-out nail.tbl \
        --ali-out nail.ali \
        ${hmm} ${fasta}
    """
}
