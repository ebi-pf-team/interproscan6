process PRINTS_RUNNER {
    label 'prints_runner'

    input:
    tuple path(fasta), val(pval), val(release), val(switches), val(postprocessing_params)

    output:
    path "prints_output.txt"
    val release
    val(postprocessing_params)

    script:
    """
    $projectDir/bin/prints/fingerPRINTScan \
    ${pval} ${fasta} ${switches} > prints_output.txt
    """

}


process PRINTS_POSTPROCESS {
    label 'analysis_parser'

    input:
    path prints_output
    val release
    val(postprocessing_params)


    output:
    path "prints_processed"
    val release

    script:
    """

    """

}


process PRINTS_PARSER {
    label 'analysis_parser'

    input:
    path "processed_prints.out"
    val release

    output:
    path "prints_parsed.json"

    script:
    """
    python3 $projectDir/scripts/prints/parser.py \
    processed_prints.out \
    ${release} \
    > prints_parsed.json
    """

}