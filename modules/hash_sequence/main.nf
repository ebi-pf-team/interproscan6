process HASH_SEQUENCE {
    input:
    val fasta

    output:
    path "parsed_sequences.json"

    script:
    """
    python3 $projectDir/scripts/hash_sequence.py -seq ${fasta} > parsed_sequences.json
    """
}
