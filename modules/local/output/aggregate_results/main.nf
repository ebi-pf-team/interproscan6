process AGGREGATE_RESULTS {
    input:
    val result_files

    output:
    path "results_aggregated.json"

    script:
    """
    mkdir -p $projectDir/results
    mkdir -p $projectDir/results/temp
    python3 $projectDir/scripts/output/aggregate_results.py "${result_files}" > results_aggregated.json
    """
}
