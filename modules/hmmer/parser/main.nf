process HMMER_PARSER {
    label 'analysis_parser'

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
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    """
}


process GENE3D_FUNFAM_PARSER {
    label 'analysis_parser'
    when:
        "${applications}".contains("${member_db}")

    input:
        path post_processed_cath_resolve_out
        val member_db
        val applications

    output:
        path "${member_db}_out.json"

    script:
    """
    python3 $projectDir/scripts/members/gene3d/assign_cath_superfamilies.py \\
        ${postprocessing_params[1]} \\
        ${postprocessing_params[2]} \\
        ${post_processed_cath_resolve_out} "${member_db}_out.json"
    """
}
