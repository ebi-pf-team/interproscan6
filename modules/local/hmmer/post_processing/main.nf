process SFLD_POST_PROCESSER {
    input:
    path out_file
    path out_dtbl
    path alignment
    val postprocessing_params  // contains [0] bin and [1] site_annotations file path
    val tsv_pro

    output:
    path "hmmer_sfld_processed"

    script:
    """
    ${postprocessing_params.last()[0]} -O ${out_file} -d ${out_dtbl} -a ${alignment} -s ${postprocessing_params.last()[1]} -o hmmer_sfld_processed
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
