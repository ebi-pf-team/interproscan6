process REVERSE_PARSE_SEQUENCE {
    input:
    tuple val(md5), val(seq_info)

    output:
    path "fasta_sequences"

    script:
    """
    python3 $projectDir/scripts/parse_sequence.py -seq_info ${seq_info} -md5 ${md5} > fasta_sequences
    """
}
