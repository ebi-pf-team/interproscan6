process PRINTS_RUNNER {
    label 'prints_runner'

    input:
    tuple path(fasta), val(pval), val(switches), val(release), val(postprocessing_params)

    output:
    path "*._.printsOutput.txt"
    val postprocessing_params

    script:
    """
    $projectDir/bin/prints/fingerPRINTScan \
    ${pval} ${fasta} ${switches} > ${release}._.printsOutput.txt
    """

}


process PRINTS_PARSER {
    label 'analysis_parser'

    input:
    path prints_output
    val postprocessing_params


    output:
    path "prints_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prints/parser.py\
    ${prints_output} \
    ${postprocessing_params} \
    > prints_parsed.json
    """

}
