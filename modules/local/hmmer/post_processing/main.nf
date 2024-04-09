process SFLD_POST_PROCESSER {
    input:
    path out_file
    path out_dtbl
    path out_alignment
    path site_annotation
    path post_process
    val tsv_pro

    output:
    path "hmmer_sfld_processed"

    script:
    """
    ${post_process} -O ${out_file} -d ${out_dtbl} -a ${out_alignment} -s ${site_annotation} > hmmer_sfld_processed
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
