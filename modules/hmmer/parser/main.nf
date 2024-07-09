process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    val postprocessing_params

    output:
    path "hmmer_parsed_${out}.json"

    script:
    """
    python3 $projectDir/scripts/hmmer/parser_out.py ${out} > hmmer_parsed_${out}.json
    """
}

process HMMER_PARSER_WITH_ALIGNMENT {
    label 'analysis_parser'

    input:
    path out
    val postprocessing_params
    path domtbl  // domtbl is needed for SFLD .c script
    val is_sfld

    output:
    path "hmmer_parsed_*"

    /*
    postprocessing_params[2] is used for parsing the domtbl
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/scripts/hmmer/${is_sfld ? "parser_out" : "parser_domtbl"}.py \\
        ${is_sfld ? "${out}" : "${domtbl}"} \\
        ${is_sfld} \\
        ${postprocessing_params[2]} \\
        > hmmer_parsed_${out}.json
    """
}

process HMMER2_PARSER {
    label 'analysis_parser'

    input:
    path out
    path fasta  // used for filtering kinase hits in SMART_FILTER_MATCHES

    output:
    path "hmmer_parsed_*"

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
