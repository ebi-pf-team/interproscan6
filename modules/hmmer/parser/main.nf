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


process GENE3D_FUNFAM_PARSER {
    label 'analysis_parser'
    when:
        "${applications}".contains("${member_db}")

    input:
        path "${release}_${hmm}*"
        path "${release}_${hmm}*"
        path "${hmm}_alignment"
        val postprocessing_params

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
