process MATCHLOOKUP {
    input:
    val fasta
    val applications

    output:
    path 'not_precalc'

    script:
    """
    python $projectDir/scripts/lookup/check_precalc.py -fasta ${fasta} -appl ${applications} > not_precalc
    """
}