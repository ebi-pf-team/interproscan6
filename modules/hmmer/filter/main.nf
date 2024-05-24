/*
These functions parse the output from the post-processing steps,
adding these data to the internal IPS6 JSON structure and also 
filtering the matches in the IPS6 JSON structure in order to 
only retain those that passed the post-processing.
*/

process PANTHER_FILTER_MATCHES {
    container 'docker.io/library/treegrafter'
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
    python3 $projectDir/scripts/members/sfld/sfld_process_post_processed.py \
        '${slfd_post_processed_output}' \
        ${ips6_json} \
        > '${ips6_json}.post.processed.json'
    """
}
