/*
These scripts implement post-processing of the HMMER files.
The output from these post-processing scripts and tools are
parsed and applied to the internal IPS6 JSON structure in the
filters module/
*/

process CATH_RESOLVE_HITS {
    label 'analysis_parser'

    input:
        path out_file
        val postprocessing_params  // [0] evalue [1] control factor
        val is_test

    output:
        path "${out_file}.cath.resolved.out"

    // cath_resolve_hits is a third party tool used to minimise suprious hits
    script:
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        touch "gene3d_out.cath.resolved.out"
        """
    else
        """
        /opt/cath-tools/cath-resolve-hits \\
            ${out_file} \\
            --input-for hmmsearch_out \\
            ${postprocessing_params[0]} > "${out_file}.cath.resolved.out"
        """
}


process FUNFAM_CATH_RESOLVE_HITS {
    /* Funfam has one hmmer.out file per cath superfamily
    Runs cath_resolve hits for each hmmer.out file as a
    single job.
    Generates an output file per hmmer.out file.
    */
    label 'analysis_parser'

    input:
        path hmmer_out_files
        val postprocessing_params  // [0] evalue [1] control factor
        val is_test

    output:
        path "*.cath.resolved.out"

    // cath_resolve_hits is a third party tool used to minimise spurious hits
    script:
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        touch "funfam.cath.resolved.out"
        """
    else
        """
        for hmmer_file in ${hmmer_out_files}
        do
            base_name=\$(basename \$hmmer_file .out)
            /opt/cath-tools/cath-resolve-hits \\
                \$hmmer_file \\
                --input-for hmmsearch_out \\
                ${postprocessing_params[0]} > "\${base_name}.cath.resolved.out"
        done
        """
}


process ADD_CATH_SUPERFAMILIES {
    label 'analysis_parser'

    input:
        path cath_resolve_out
        val postprocessing_params
        val is_test
    /*
    Post-processing params:
    1. cath_resolve_hits_switches
    2. model2sf_map - path to data file
    3. discontinuous_regs - path to data file
    4. assign_cath_superfamilies - path to py script
    */

    output:
        path "${cath_resolve_out}.cath_superfamilies"

    script:
    if ( is_test )
        """
        echo ${postprocessing_params[1]}
        echo ${postprocessing_params[2]}
        echo ${postprocessing_params[3]}
        touch "${cath_resolve_out}.cath_superfamilies"
        """
    else
        """
        python3 ${postprocessing_params[3]} \\
            ${postprocessing_params[1]} \\
            ${postprocessing_params[2]} \\
            ${cath_resolve_out} "${cath_resolve_out}.cath_superfamilies"
        """
}


process HAMAP_POST_PROCESSER {
    label 'analysis_parser'

    input:
        val postprocessing_params
        path fasta
        path tlb
        val is_test
    /*
    post-processing params:
    0. models dir
    1. flags for pfsearch
    */

    output:
        path "hamap_pfsearch_output"

    /*
    Args for pfsearch_wrapper.py
    hmmer_tbl_path -- hmmer tlb file
    fasta_file -- input fasta
    fasta_filtered_file -- output fasta file
    output_file -- another output file
    model_dir -- path to dir containing hamap profiles for pfsearch
    path to pfsearchV3 bin
    flags for pfsearch
    */
    script:
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        echo ${postprocessing_params[1]}
        touch "hamap_pfsearch_output"
        """
    else
        """
        python3 $projectDir/interproscan/scripts/members/hamap/pfsearch_wrapper.py \
            ${tlb} \
            ${fasta} \
            "seqs_with_hits.faa" \
            "hamap_pfsearch_output" \
            ${postprocessing_params[0]} \
            "/opt/pftools/pfsearchV3" \
            ${postprocessing_params[1]} > "print.statements.out"
        """
}


process PANTHER_POST_PROCESSER {
    label 'treegrafter_analysis'

    input:
        path out_file
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        path fasta
        val is_test

    output:
        path "treegrafter_processed_panther_hits"

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
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        echo ${postprocessing_params[1]}
        touch treegrafter_processed_panther_hits
        """
    else
        """
        mkdir data
        python3 $projectDir/interproscan/scripts/members/panther/treegrafter.py \\
            run \\
            ${fasta} \\
            ${out_file} \\
            ${postprocessing_params[0]} \\
            -e ${postprocessing_params[1]} \\
            -o treegrafter_processed_panther_hits \\
            --epa-ng /opt/epa-ng/bin/epa-ng \\
            --keep
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
        path out_file
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        val member
        path alignment
        path out_dtbl
        val is_test

    output:
        path "${out_file}.processed.out"

    script:
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        touch "${out_file}.processed.out"
        """
    else
        """
        $projectDir/interproscan/bin/sfld/sfld_postprocess \
            --alignments '${alignment}' \
            --dom '${out_dtbl}' \
            --hmmer-out '${out_file}' \
            --site-info '${postprocessing_params[0]}' \
            --output '${out_file}.processed.out'
        """
}


process SUPERFAMILY_POST_PROCESSER {
    label 'analysis_parser'

    input:
        path hmmscan_out
        val postprocessing_params
        path fasta
        val is_test

    /*
    postprocessing_params[0] = bin
    postprocessing_params[1] = self_hits
    postprocessing_params[2] = cla
    postprocessing_params[3] = model_tab
    postprocessing_params[4] = pdbj95d
    postprocessing_params[5] = binary_switches
    */

    output:
        path "${hmmscan_out}_ass3_output"

    script:
    if ( is_test )
        """
        echo ${postprocessing_params[0]}
        echo ${postprocessing_params[1]}
        echo ${postprocessing_params[2]}
        echo ${postprocessing_params[3]}
        echo ${postprocessing_params[4]}
        echo ${postprocessing_params[5]}
        touch "${hmmscan_out}_ass3_output"
        """
    else
        """
        perl ${postprocessing_params[0]} \\
        -s ${postprocessing_params[1]} \\
        -r ${postprocessing_params[2]} \\
        -m ${postprocessing_params[3]} \\
        -p ${postprocessing_params[4]} \\
        ${postprocessing_params[5]} \\
        ${fasta} ${hmmscan_out} ${hmmscan_out}_ass3_output
        """
}
