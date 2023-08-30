process WRITERESULTS {
    publishDir "${params.projectDir}/results", mode: 'copy'

    input:
    val xref_results
    val formats
    val output_path

    script:
    """
    python $projectDir/scripts/write_output.py --results ${xref_results} --format ${formats} --output_path ${output_path}
    """
}