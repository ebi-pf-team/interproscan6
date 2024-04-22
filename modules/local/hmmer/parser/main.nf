process HMMER_PARSER {
    input:
    path out
    path domtbl
    path alignment  // Needed for SFLD post processing
    val postprocessing_params
    val tsv_pro
    val sites  // SFLD and CDD predict sites

    output:
    path "hmmer_parsed_*"

    script:
    """
    rm -f ${alignment}
    python3 $projectDir/scripts/hmmer/${tsv_pro ? "parser_out" : "parser_domtbl"}.py \\
    ${tsv_pro ? "${out}" : "${domtbl}"} \\
    ${sites ? "true" : "false"} \\
    > hmmer_parsed_${out}.json
    """
}
