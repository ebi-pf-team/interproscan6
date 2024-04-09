process SFLD_POST_PROCESSER {
    input:
    path out
    val tsv_pro

    output:
    path "hmmer_sfld_processed"

    script:
    """
    ... > ???
    """
}


process FUNFAM_PSOT_PROCESSER {
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
