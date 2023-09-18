process WRITERESULTS {
    publishDir "${params.projectDir}/results", mode: 'copy'

    input:
    val collected_outputs
    val format
    val output_path

    script:
    """
    cat ${collected_outputs.join(" ")} > $projectDir/${output_path}.tmp
    python3 $projectDir/scripts/write_output.py --results $projectDir/${output_path}.tmp --format ${format} --output_path $projectDir/${output_path}    """
}