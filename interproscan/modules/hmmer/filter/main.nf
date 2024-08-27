/*
These functions parse the output from the post-processing steps,
adding these data to the internal IPS6 JSON structure and also
filtering the matches in the IPS6 JSON structure in order to
only retain those that passed the post-processing.
*/

process FUNFAM_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_jsons
        path cath_resolved_outs
        val postprocessing_params
    // postprocessing_params[6] funfam release

    output:
        path "*.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/funfam/run_match_filtering.py \\
        ${ips6_jsons} \\
        ${cath_resolved_outs} \\
        ${postprocessing_params[6]} > debug
    """
}


process GENE3D_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path cath_resolve_out_with_superfams
        path ips6_json
        val postprocessing_params

    output:
        path "${ips6_json}.post.processed.json"
        path "cath.superfamilies"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/gene3d/filter_ips6_hits.py \\
        ${ips6_json} \\
        ${cath_resolve_out_with_superfams} \\
        ${ips6_json}.post.processed.json \\
        cath.superfamilies \\
        ${postprocessing_params[4]}
    """
}


process HAMAP_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path pfsearch_wrapper_output

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/hamap/filter_ips6_hits.py \\
        ${ips6_json} \\
        ${pfsearch_wrapper_output} \\
        ${ips6_json}.post.processed.json
    """
}


process PANTHER_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path treegrafter_output
        val postprocessing_params

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/panther/process_treegrafter_hits.py \\
        ${treegrafter_output} \\
        ${ips6_json} \\
        ${postprocessing_params[2]} \\
        ${ips6_json}.post.processed.json
    """
}


process PFAM_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        val postprocessing_params

    output:
        path "${ips6_json}.post.processed.json"

    /*
    CL input to postprocess.py:
    0. Str repr of the path to the internal IPS6 JSON of raw PFAM results
    1. The min length specified in post-processing params in the members.config file
    2. Str repr of the path to the PFAM seed file
    3. Str repr of the path to the PFAM clan file
    4. Str repr of the path to the PFAM dat file
    5. Str repr for the output file path
    */
    script:
    """
    python3 $projectDir/interproscan/scripts/members/pfam/postprocess.py \\
        ${ips6_json} \\
        '${postprocessing_params[0]}' \\
        '${postprocessing_params[1]}' \\
        '${postprocessing_params[2]}' \\
        '${postprocessing_params[3]}' \\
        '${ips6_json}.post.processed.json'
    """
}


process PIRSF_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path dtbl_file  // needed to get the sequence length
        val postprocessing_params   // [0] path to the PIRSF.dat file

    /* PIRSF uses the total sequence len, tlen in the hmmscan dtbl output,
    for the HmmLength */
    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/pirsf/filter_ips6_hits.py \
        ${ips6_json} \
        ${dtbl_file} \
        '${postprocessing_params[0]}' \
        '${ips6_json}.post.processed.json'
    """
}


process PIRSR_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        val postprocessing_params   // [0] path to rules file

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/pirsr/filter.py \\
        ${ips6_json} \\
        ${postprocessing_params[0]} \\
        ${ips6_json}.post.processed.json
    """
}


process SFLD_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path slfd_post_processed_output

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/sfld/sfld_process_post_processed.py \\
        ${slfd_post_processed_output} \\
        ${ips6_json} > ${ips6_json}.post.processed.json
    """
}


process SMART_FILTER_MATCHES {
    label 'analysis_parser'
    /*
    It needs the FASTA file becauce when both Ser-Thr and Tyr
    kinase matches are found in a sequence, the domains
    are checked again using regex checks against the protein sequence.
    */

    input:
        path ips6_json
        path fasta

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/smart/filter_ips6_hits.py \\
        ${ips6_json} \\
        ${fasta} \\
        ${ips6_json}.post.processed.json
    """
}


process SUPERFAMILY_FILTER_MATCHES {
    // parses the output from SUPERFAMILY_POST_PROCESSER
    label 'analysis_parser'

    input:
    path ass3_out
    path hmm_lib

    output:
    path "superfamily_parsed_*"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/superfamily/parse_superfamily_out.py \\
        ${hmm_lib} \\
        ${ass3_out} \\
        superfamily_parsed_${ass3_out}.json
    """
}
