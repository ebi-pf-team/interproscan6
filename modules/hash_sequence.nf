process XREFS {
    input:
    val fasta

    output:
    path 'parsed_sequences.json'

    script:
    """
    python3 $projectDir/scripts/sequences_parse.py -fasta ${fasta} > parsed_sequences.json
    """
}