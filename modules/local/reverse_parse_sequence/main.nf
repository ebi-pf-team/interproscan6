process REVERSE_PARSE_SEQUENCE {
    input:
    val input_file
    val seq_info

    output:
    path "fasta_sequences", optional: true

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py -file ${input_file} --seq_info ${seq_info} > fasta_sequences
    """
}
