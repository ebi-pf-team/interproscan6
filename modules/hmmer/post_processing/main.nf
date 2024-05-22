
process PANTHER_POST_PROCESSER {
    container 'docker.io/library/treegrafter'
    label 'analysis_parser'

    input:
        path ips6_json
        path out_file
        path out_dtbl
        val postprocessing_params // [0] evalue [1] control factor
        val tsv_pro
        path fasta

    output:
        path "${out_dtbl}.post.processed.json"

    /*
    Input args for TreeGrafter:
    1. input FASTA file
    2. Out file from HMMSearch
    3. TreeGrafter data dir (from members config)
    4. -e is the evalue threshold (from members.config)
    5. -o output file to be created
    6. --keep temp directory
    */
    script:
    """
    mkdir data
    python3 $projectDir/scripts/members/panther/treegrafter.py \
        run \
        ${fasta} \
        ${out_file} \
        ${postprocessing_params[0]} \
        -e ${postprocessing_params[1]} \
        -o processed_panther_hits \
        --epa-ng /epa-ng/bin/epa-ng \
        --keep

    python3 $projectDir/scripts/members/panther/process_treegrafter_hits.py \
        processed_panther_hits \
        "${out_dtbl}.post.processed.json"
    """
}


process SFLD_POST_PROCESSER {
    /*
    Runs an in-house post-processing C script that filters the SFLD hits to
    retain only those where all sites are conserved, and add site data to
    the remaining hits. The out file is then parsed and the results in the IPS6
    json file are filtered and site data is added to the remaining hits.
    */
    label 'analysis_parser'

    input:
        path ips6_json
        path out_file
        path out_dtbl
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        path alignment
        val tsv_pro

    output:
        path "${tsv_pro ? "${out_file}.post.processed.out.json" : "${out_dtbl}.post.processed.dtbl.json"}"

    script:
        """
        ${postprocessing_params[0]} \
            --alignments '${alignment}' \
            --dom '${out_dtbl}' \
            --hmmer-out '${out_file}' \
            --site-info '${postprocessing_params[1]}' \
            --output '${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}'

        python3 $projectDir/scripts/members/sfld/sfld_process_post_processed.py \
            '${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}' \
            ${ips6_json} \
            > '${tsv_pro ? "${out_file}.post.processed.out.json" : "${out_dtbl}.post.processed.dtbl.json"}'
        """
}
