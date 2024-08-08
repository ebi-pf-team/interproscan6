process PARSE_SEQUENCE {
    label 'io'

    input:
    val fasta_file
    val applications

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/interproscan/scripts/parse_sequence.py ${fasta_file} ${applications} > parsed_sequences
    """
}