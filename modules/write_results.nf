process WRITERESULTS {
    input:
    val collected_sequences
    val collected_outputs
    val format
    val output_path

    script:
    """
    cat ${collected_sequences.join(" ")} > $projectDir/results/temp/sequences_hash.tmp
    cat ${collected_outputs.join(" ")} > $projectDir/results/temp/matches_result.tmp
    python3 $projectDir/scripts/write_output.py -matches $projectDir/results/temp/matches_result.tmp -seq $projectDir/results/temp/sequences_hash.tmp -format ${format} -out $projectDir/results/${output_path}
    """
}

