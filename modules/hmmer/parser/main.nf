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
    path out
    path domtbl
    val postprocessing_params
    path(alignment), optional: true

    /*
    member_db --> true/false is to tell the domtbl parser if to retrieve site data
    postprocessing_params[2] is used for panther when parsing the domtbl 
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/scripts/hmmer/${tsv_pro ? "parser_out" : "parser_domtbl"}.py \\
        ${tsv_pro ? "${out}" : "${domtbl}"} \\
        ${member_db == "sfld" ? "true" : "false"} \\
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    if [${member_db} != "sfld" || ${member_db} != "gene3d" || ${member_db} != "funfam"]; then
        rm ${alignment}
    fi
    """
}
