process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    path domtbl
    val postprocessing_params
    val tsv_pro
    val get_sites

    output:
    path "hmmer_parsed_*"

    /*
    get_sites --> "true" or "false" is to tell the domtbl parser if to retrieve site data
    postprocessing_params[2] is used for panther when parsing the domtbl 
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/scripts/hmmer/${tsv_pro ? "parser_out" : "parser_domtbl"}.py \\
        ${tsv_pro ? "${out}" : "${domtbl}"} \\
        ${get_sites} \\
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    """
}


process HMMER2_PARSER {
    label 'analysis_parser'

    input:
    path out
    val postprocessing_params

    output:
    path "hmmer_parsed_*"
    val postprocessing_params

    /*
    get_sites --> "true" or "false" is to tell the domtbl parser if to retrieve site data
    postprocessing_params[2] is used for panther when parsing the domtbl 
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/scripts/hmmer/parse_hmmpfam_out.py \\
        ${out} \\
        > hmmer_parsed_${out}.json
    """
}
