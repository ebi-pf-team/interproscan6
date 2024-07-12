process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out

    output:
    path "hmmer_parsed_${out}.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/hmmer/parser_out.py ${out} > hmmer_parsed_${out}.json
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
    python3 $projectDir/interproscan/scripts/hmmer/parse_hmmpfam_out.py \\
        ${out} \\
        > hmmer_parsed_${out}.json
    """
}
