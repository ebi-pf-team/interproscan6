process UNION_RESULTS {
    input:
    val pre_calc_results
    val analysis_results

    output:
    path "all_results.json"

    script:
    """
    python3 $projectDir/scripts/post_proc/union_results.py -pre_calc "${pre_calc_results}" -analysis "${analysis_results}" > all_results.json
    """
}
