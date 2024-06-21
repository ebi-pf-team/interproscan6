process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    val postprocessing_params // used to post processing

    output:
    path "hmmer_parsed_*"
    val postprocessing_params

    script:
    """
    python3 $projectDir/scripts/hmmer/parser_out.py ${out} > hmmer_parsed_${out}.json
    """
}

process HMMER_PARSER_WITH_ALIGNMENT {
    label 'analysis_parser'

    input:
        path out
        path domtbl  // domtbl is needed for SFLD .c script
        path alignment
        val postprocessing_params


        output:
        path "hmmer_parsed_*"
        path domtbl
        val postprocessing_params
        path(alignment)
        val is_sfld  // just SFLD needs to retrieve site data and parse dtbl
        /*
        postprocessing_params[2] is used for panther when parsing the domtbl
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


process HMMER_PARSER_TBL {
    label 'analysis_parser'

    input:
        path out
        path tbl  // used to post processing
        val postprocessing_params // used to post processing

        output:
        path "hmmer_parsed_*"
        path tbl
        val postprocessing_params

        script:
        """
        python3 $projectDir/scripts/hmmer/parser_out.py ${out} > hmmer_parsed_${out}.json
        """
}
