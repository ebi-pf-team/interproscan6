process CDD {
    input:
    val fasta
    val applications

    output:
    path 'parsed_cdd'

    script:
    """
    python $projectDir/scripts/homologs_search/cdd.py -fasta ${fasta} -appl ${applications} > parsed_cdd
    """
}