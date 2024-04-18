process GET_ORFS {
    container 'docker.io/library/easel'

    input:
    val fasta_file
    val(strand)
    val(methionine)
    val(min_len)
    val(genetic_code)

    output:
    path "orfs_sequences.fasta"

    script:
    """
    /opt/easel/easel/miniapps/esl-translate \
        -c ${genetic_code} \
        -l ${min_len} \
        ${fasta_file} \
        > orfs_sequences.fasta
    """
}
