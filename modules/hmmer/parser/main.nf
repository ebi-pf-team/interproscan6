process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    path domtbl
    val postprocessing_params
    val get_sites
    val member_db

    output:
    path "hmmer_parsed_*"

    /*
    get_sites --> "true" or "false" is to tell the domtbl parser if to retrieve site data
    postprocessing_params[2] is used for panther when parsing the domtbl
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/scripts/hmmer/parser_out.py \\
        ${out} \\
        ${get_sites} \\
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    """
}
