process WRITE_RESULTS {
    label 'io'
    publishDir "$projectDir/results/", mode: 'copy'

    input:
    val sequences
    path matches
    val format
    val output_path

    script:
    """
    cat ${sequences.join(" ")} > $projectDir/results/temp/sequences_hash.tmp
    python3 $projectDir/scripts/output/write_results.py $projectDir/results/temp/sequences_hash.tmp ${matches} ${format} $projectDir/${output_path} > debug_out
    """
}
