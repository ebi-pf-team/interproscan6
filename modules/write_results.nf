process WRITERESULTS {
    publishDir "${params.projectDir}/results", mode: 'copy'

    input:
    val xref_results
    val format
    val output_path

    script:
    """
    python3 $projectDir/scripts/write_output.py --results ${xref_results} --format ${format} --output_path $projectDir/${output_path}
    """
}