process PRINTS_RUNNER {
    label 'prints_runner'

    input:
    tuple path(fasta), val(pval), val(switches), val(release), val(postprocessing_params)

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
    path "prints_processed.out"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prints/postprocess.py\
    prints_output.txt\
    ${postprocessing_params[0]}\
    ${release} \
    > prints_processed.out
    """

}


process PRINTS_PARSER {
    label 'analysis_parser'

    input:
    path "processed_prints.out"

    output:
    path "prints_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prints/parser.py \
    processed_prints.out \
    > prints_parsed.json
    """

}