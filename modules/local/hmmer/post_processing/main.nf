process SFLD_POST_PROCESSER {
    input:
    path out_file
    path out_dtbl
    path alignment
    val postprocessing_params  // contains [0] bin and [1] site_annotations file path
    val tsv_pro

    output:
    path touch "${out_file}.processed.out"
    path touch "${out_dtbl}.processed.dtbl"
    path alignment
    val postprocessing_params

    script:
    """
    ${postprocessing_params[0]} \\
        ${tsv_pro ? -O ${out_file} : -d ${out_dtbl}}
        -a ${alignment} \\
        -s ${postprocessing_params[1]} \\
        -o ${tsv_pro ? -O "${out_file}.processed.out" : -d "${out_dtbl}.processed.dtbl"}
    
    if [tsv_pro]; then
        touch "${out_dtbl}.processed.dtbl"
    else
        touch "${out_file}.processed.out"
    fi
    """
}


process FUNFAM_POST_PROCESSER {
    input:
    path out
    val tsv_pro

    output:
    path "hmmer_funfam_processed"

    script:
    """
    ... > ???
    """
}


process GENE3D_POST_PROCESSER {
    input:
    path out
    val tsv_pro

    output:
    path "hmmer_gene3d_processed"

    script:
    """
    ... > ???
    """
}
