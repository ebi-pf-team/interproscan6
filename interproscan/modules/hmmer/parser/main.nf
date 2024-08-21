/*
    This processes invoke scripts that parser the output from HMMER
    (hmmsearch, hmmscan, hmmpfam) into the IPS6 JSON structure.
*/
process HMMER_PARSER {
    label 'analysis_parser'

    input:
    path out
    val member

    output:
    path "hmmer_parsed_${out}.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/hmmer/parser_out.py \\
        ${out} \\
        ${member} \\
        hmmer_parsed_${out}.json
    """
}


process FUNFAM_HMMER_PARSER {
    /*
    Parses output from HMMER3, but takes in multiple hmmer.out files, and generates
    an JSON file for each --> output all generated json files
    */
    label 'analysis_parser'

    input:
    path hmmer_out_files
    val member

    output:
    path "hmmer_parsed_*.json"

    script:
    """
    for hmmer_file in ${hmmer_out_files}
    do
        base_name=\$(basename \$hmmer_file .out)
        python3 $projectDir/interproscan/scripts/hmmer/parser_out.py \$hmmer_file ${member} hmmer_parsed_\${base_name}.json
    done
    """
}


process HMMER2_PARSER {
    label 'analysis_parser'

    input:
    path out
    path fasta  // used for filtering kinase hits in SMART_FILTER_MATCHES
    val member

    output:
    path "hmmer_parsed_${out}.json"

    /*
    get_sites --> "true" or "false" is to tell the domtbl parser if to retrieve site data
    postprocessing_params[2] is used for panther when parsing the domtbl
    These won't be needed in the python script call when using only HMMER.out
    */
    script:
    """
    python3 $projectDir/interproscan/scripts/hmmer/parse_hmmpfam_out.py \\
        ${out} \\
        ${member} \\
        hmmer_parsed_${out}.json
    """
}


process HMMER_SCAN_PARSER {
    label 'analysis_parser'

    input:
    path out
    val member

    output:
    path "hmmer_parsed_${out}.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/hmmer/parser_scan_out.py \\
        ${out} \\
        ${member} \\
        hmmer_parsed_${out}.json
    """
}
