process AGGREGATE_RESULTS {
    label 'io'
    
    input:
    val result_files  //  paths to the json files created in SEQUENCE_ANALYSIS concatenated with parsed_matches (lookup)

    output:
    path "results_aggregated.json"

    script:
    """
    mkdir -p $projectDir/results
    mkdir -p $projectDir/results/temp
    python3 $projectDir/scripts/output/aggregate_results.py "${result_files}" > results_aggregated.json
    """
}
