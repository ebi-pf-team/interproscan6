process WRITE_RESULTS {
    publishDir "$projectDir/results/", mode: 'copy'

    input:
    val collected_sequences
    path matches
    val format
    val output_path

    script:
    """
    cat ${collected_sequences.join(" ")} > $projectDir/results/temp/sequences_hash.tmp
    python3 $projectDir/scripts/write_output.py -matches ${matches} -seq $projectDir/results/temp/sequences_hash.tmp -format ${format} -out $projectDir/results/${output_path}
    """
}

