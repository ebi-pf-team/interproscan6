process PARSE_SEQUENCE {
    input:
    val fasta_file

    output:
    path "parsed_sequences"

    script:
    """
    if [[ ! -f ${fasta_file} ]]; then
        echo "File not found: ${fasta_file}"
        exit 1
    fi
    python3 $projectDir/scripts/parse_sequence.py ${fasta_file} > parsed_sequences
    """
}
