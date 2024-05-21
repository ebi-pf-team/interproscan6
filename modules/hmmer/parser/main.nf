process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    path domtbl
    path alignment  // Needed for SFLD post processing
    val postprocessing_params
    val tsv_pro
    val member_db

    output:
    path "hmmer_parsed_*"
    val postprocessing_params
    path(alignment), optional: true

    script:
    """
    python3 $projectDir/scripts/hmmer/${tsv_pro ? "parser_out" : "parser_domtbl"}.py \\
        ${tsv_pro ? "${out}" : "${domtbl}"} \\
        ${sites ? "true" : "false"} \\
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    if [${member_db} == "antifam" || ${member_db} == "ncbifam"]; then
        rm ${alignment}
    fi
    """
}
