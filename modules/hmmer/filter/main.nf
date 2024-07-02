/*
These functions parse the output from the post-processing steps,
adding these data to the internal IPS6 JSON structure and also
filtering the matches in the IPS6 JSON structure in order to
only retain those that passed the post-processing.
*/

process FUNFAM_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path cath_resolved_out
        val postprocessing_params
    /*
    Post-processing params are needed when cath_resolved hits
    process feeds into add_cath_superfamilies process
    */
    
    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/scripts/members/funfam/filter_ips6_hits.py \\
        ${ips6_json} \\
        ${cath_resolved_out} \\
        ${ips6_json}.post.processed.json \\
        ${postprocessing_params[6]}
    """
}


process GENE3D_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path cath_resolve_out_with_superfams
        path ips6_json
    
    output:
        path "${ips6_json}.post.processed.json"
        path "cath.superfamilies"

    script:
    """
    python3 $projectDir/scripts/members/gene3d/filter_ips6_hits.py \\
        ${ips6_json} \\
        ${cath_resolve_out_with_superfams} \\
        ${ips6_json}.post.processed.json \\
        cath.superfamilies
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
    python3 $projectDir/scripts/members/hamap/filter_ips6_hits.py \\
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
    python3 $projectDir/scripts/members/panther/process_treegrafter_hits.py \
        ${treegrafter_output} \
        ${ips6_json} \
        ${postprocessing_params[2]} > ${ips6_json}.post.processed.json
    """
}

process SFLD_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path slfd_post_processed_output
        val postprocessing_params

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/scripts/members/sfld/sfld_process_post_processed.py \\
        ${slfd_post_processed_output} \\
        ${ips6_json} > ${ips6_json}.post.processed.json
    """
}

process PFAM_FILTER_MATCHES {
    label 'analysis_parser'

    input:
        path ips6_json
        path out_file
        path out_dtbl
        val postprocessing_params
        path alignment

    output:
        path "${ips6_json}.post.processed.json"

    script:
    """
    python3 $projectDir/scripts/members/pfam/postprocess.py \
        --hmm_parsed ${ips6_json} \
        --min_length '${postprocessing_params[0]}' \
        --seed '${postprocessing_params[1]}' \
        --clans '${postprocessing_params[2]}' \
        --dat '${postprocessing_params[3]}' \
        > '${ips6_json}.post.processed.json'
    """
}
