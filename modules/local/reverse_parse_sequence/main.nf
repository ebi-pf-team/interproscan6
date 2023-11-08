process REVERSE_PARSE_SEQUENCE {
    input:
    val md5
    val seq_info

    output:
    path "fasta_sequences"

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py -md5 ${md5} -seq_info ${seq_info} > fasta_sequences
    """
}
