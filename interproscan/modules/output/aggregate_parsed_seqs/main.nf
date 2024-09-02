process AGGREGATE_PARSED_SEQS {
    label 'io'
    
    input:
    val parsed_sequences_files

    output:
    path "parsed_seqs_aggregated.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/output/aggregate_parsed_seqs.py \\
        "${parsed_sequences_files}" \\
        parsed_seqs_aggregated.json
    """
}
