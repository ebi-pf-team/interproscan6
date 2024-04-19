process WRITE_RESULTS {
    publishDir "$projectDir/results/", mode: 'copy'

    input:
    val sequences
    path matches
    val format
    val output_path

    script:
    """
    python3 $projectDir/scripts/output/write_results.py ${sequences} ${matches} ${format} $projectDir/${output_path} > debug_out
    """
}
