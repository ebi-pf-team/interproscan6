/*
These scripts implement post-processing of the HMMER files.
The output from these post-processing scripts and tools are 
parsed and applied to the internal IPS6 JSON structure in the 
filters module/
*/

process CATH_RESEOLVE_HITS {
    label 'analysis_parser'

    input:
        path out_file
        val postprocessing_params  // [0] evalue [1] control factor

    output:
        path "${out_file}.cath.resolved.out"
        val postprocessing_params

    // cath_resolve_hits is a third party tool used to minimise suprious hits
    script:
    """
    /opt/cath-tools/cath-resolve-hits \\
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
        path fasta
        path tlb
        val postprocessing_params
    /*
    post-processing params:
    0. models dir
    1. flags for pfserach
    */

    output:
        path "hamap_pfsearch_output"

    /*
    hmmer_tbl_path -- hmmer tlb file
    fasta_file -- input fasta
    fasta_filtered_file -- output fasta file
    output_file -- another output file
    model_dir -- path to dir containing hamap profiles for pfsearch
    */
    script:
    """
    python3 $projectDir/scripts/members/hamap/pfsearch_wrapper.py \
        ${tlb} \
        ${fasta} \
        "seqs_with_hits.faa" \
        "hamap_pfsearch_output" \
        ${postprocessing_params[0]} \
        "/opt/pftools/var/lib/pftools/bin/pfsearchV3" \
        ${postprocessing_params[1]} > "print.statements.out"
    """
}


process PANTHER_POST_PROCESSER {
    label 'analysis_parser'

    input:
        path out_file
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        path fasta

    output:
        path "treegrafter_processed_panther_hits"
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
        -o treegrafter_processed_panther_hits \
        --epa-ng /opt/epa-ng/bin/epa-ng \
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
        path out_dtbl
        val postprocessing_params // contains [0] bin and [1] site_annotations file path
        path alignment
        val tsv_pro

    output:
        path "${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}"
        val postprocessing_params

    script:
        """
        ${postprocessing_params[0]} \
            --alignments '${alignment}' \
            --dom '${out_dtbl}' \
            --hmmer-out '${out_file}' \
            --site-info '${postprocessing_params[1]}' \
            --output '${tsv_pro ? "${out_file}.processed.out" : "${out_dtbl}.processed.dtbl"}'
        """
}
