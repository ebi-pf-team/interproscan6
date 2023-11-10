process PARSE_SEQUENCE {
    input:
    val fasta_file

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py -fasta ${fasta_file} > parsed_sequences
    """
}
