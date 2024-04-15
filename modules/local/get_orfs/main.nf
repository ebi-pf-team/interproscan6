process GET_ORFS {
    input:
    val fasta_file

    output:
    path "orfs_sequences.fasta"

    script:
    """
    ${params.orfs.translate} -i ${fasta_file} -o orfs_sequences.fasta
    """
}
