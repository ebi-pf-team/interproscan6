process AGGREGATE_RESULTS {
    label 'io'
    
    input:
    val result_files  //  paths to the json files created in SEQUENCE_ANALYSIS concatenated with parsed_matches (lookup)

    output:
    path "results_aggregated.json"

    /* results_files can become longer than the bash limit resulting in 
    a 126 exit status and 'Arugment list too long' error. Therefore, the
    str respresentation of a list that is results_files is 
    broken up into chunks, and each chunk is added to the final 
    results_aggregated.json file. */
    script:
    """
    mkdir -p $projectDir/results
    mkdir -p $projectDir/results/temp
    echo "{}" > results_aggregated.json
    IFS=',' read -ra paths <<< "${result_files}"
    for ((i = 0; i < ${#paths[@]}; i += 20)); do
        python3 $projectDir/interproscan/scripts/output/aggregate_results.py \\
        "\${batch[@]}" \\
        results_aggregated.json
    """
}
