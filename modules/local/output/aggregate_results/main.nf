process AGGREGATE_RESULTS {
    input:
    val result_files

    output:
    path "results_aggregated.json"

    script:
    """
    python3 $projectDir/scripts/output/aggregate_results.py "${result_files}" > results_aggregated.json
    """
}
