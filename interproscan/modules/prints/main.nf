process PRINTS_RUNNER {
    label 'prints_runner'

    input:
    tuple path(fasta), val(hierarchy), val(pval), val(switches)

    output:
    path "*._.printsOutput.txt"
    val hierarchy

    script:
    """
    $projectDir/bin/prints/fingerPRINTScan \
    ${pval} ${fasta} ${switches} > printsOutput.txt
    """

}


process PRINTS_PARSER {
    label 'analysis_parser'

    input:
    path prints_output
    val hierarchy


    output:
    path "prints_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prints/parser.py \\
    ${prints_output} \\
    ${hierarchy} \\
    prints_parsed.json
    """

}
