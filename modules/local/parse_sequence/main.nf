process PARSE_SEQUENCE {
    label 'io'

    input:
    val fasta_file

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py ${fasta_file} > parsed_sequences
    """
}
