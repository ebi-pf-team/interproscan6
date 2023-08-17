process WRITERESULTS {
    publishDir "${params.projectDir}/results", mode: 'copy'

    input:
    val result_parsed
    val formats

    script:
    """
    python $projectDir/scripts/write_output.py --results ${result_parsed} --format ${formats}
    """
}