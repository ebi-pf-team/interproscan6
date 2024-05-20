
process CATH_RESEOLVE_HITS {
    container 'docker.io/library/cathtools'
    label 'analysis_parser'

    input:
        path out_file
        path out_dtbl
        path alignment  // note used here. Needed for the SFLD pipeline
        val postprocessing_params  // [0] evalue [1] control factor

    output:
        path "${out_file}.cath.resolved.out"
        val postprocessing_params

    // cath_resolve_hits is a third party tool used to minimise suprious hits
    script:
    """
    rm -f ${alignment}
    /cath-tools/cath-resolve-hits \\
        ${out_file} \\
        --input-for hmmsearch_out \\
        ${postprocessing_params[0]} > "${out_file}.cath.resolved.out"
    """
}


process ADD_CATH_SUPERFAMILIES {
    label 'analysis_parser'

    input:
        path cath_resolve_out
        val postprocessing_params
        val member_db

    output:
        path "${out_file}.${member_db}.json"
        val postprocessing_params

    script:
    """
    python3 $projectDir/members/gene3d_funfam/assign_cath_superfamilies.py \\
        ${postprocessing_params[1]} \\
        ${postprocessing_params[2]} \\
        ${cath_resolve_out} "${out_file}.${member_db}.json"
    """
}


process PANTHER_POST_PROCESSER {
    container 'docker.io/library/treegrafter'
    label 'analysis_parser'

    input:
        path out_file
        path out_dtbl
        path alignment  // note used here. Needed for the SFLD pipeline
        val postprocessing_params  // [0] evalue [1] control factor
        path fasta

    output:
        path out_file
        path "${out_dtbl}.post_processed.dtbl"
        path alignment
        val postprocessing_params

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
        ${out_dtbl}
    """
}


process SFLD_POST_PROCESSER {
    label 'analysis_parser'

    input:
        path out_file
        path out_dtbl
        path alignment
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        val tsv_pro

    output:
        path "${out_file}.post_processed.out"
        path "${out_dtbl}.post_processed.dtbl"
        path alignment
        val postprocessing_params // need to make up input var nums in hmmer_parser module

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
            ${tsv_pro ? "-O '${out_file}'" : "-d '${out_dtbl}'"}

        touch ${out_file}.post_processed.out
        touch ${out_dtbl}.post_processed.dtbl
        """
}
