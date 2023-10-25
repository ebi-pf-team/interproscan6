process PARSE_SEQUENCE {
    input:
    val input_file
    val reverse

    output:
    path "parsed_sequences"

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py -file ${input_file} -reverse ${reverse} > parsed_sequences
    """
}
