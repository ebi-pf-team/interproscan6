process AGGREGATE_RESULTS {
    input:
    val result_files  //  paths to the json files created in SEQUENCE_ANALYSIS

    output:
    path "results_aggregated.json"

    script:
    """
    mkdir -p $projectDir/results
    mkdir -p $projectDir/results/temp
    python3 $projectDir/scripts/output/aggregate_results.py "${result_files}" > results_aggregated.json
    """
}
