process WRITE_RESULTS {
    publishDir "$projectDir/results/", mode: 'copy'

    input:
    val sequences
    path matches
    val format
    val output_path
    val version

    script:
    """
    cat ${sequences.join(" ")} > $projectDir/results/temp/sequences_hash.tmp
    python3 $projectDir/scripts/output/write_results.py $projectDir/results/temp/sequences_hash.tmp ${matches} ${format} $projectDir/${output_path} $version > debug_out
    """
}
