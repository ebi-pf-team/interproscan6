process GET_ORFS {
    container 'docker.io/library/easel'

    input:
    val fasta_file
    val strand
    val methionine
    val min_len
    val genetic_code

    output:
    path "orfs_sequences.fasta"

    script:
    def analysed_strand = (strand == "both") ? "" : (strand == "plus") ? "--watson" : "--crick"

    """
    echo ${strand}
    echo ${analysed_strand}
    /opt/easel/easel/miniapps/esl-translate \\
    -c ${genetic_code} \\
    -l ${min_len} \\
    ${analysed_strand} \\
    ${fasta_file} \\
    > orfs_sequences.fasta
    """
}
