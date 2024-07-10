process GET_ORFS {
    label 'get_orfs'

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
    /opt/hmmer-3.3/easel/miniapps/esl-translate \\
    -c ${genetic_code} \\
    -l ${min_len} \\
    ${analysed_strand} \\
    ${fasta_file} \\
    > orfs_sequences.fasta
    """
}
