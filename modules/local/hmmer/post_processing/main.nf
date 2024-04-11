process SFLD_POST_PROCESSER {
    input:
    path out_file
    path out_dtbl
    val sfld_params  // contains bin and site_annotations file path
    val tsv_pro

    output:
    path "hmmer_sfld_processed"

    script:
    """
    if [ -f "hmmer_alignment" ]; then
        ${sfld_params.last()[0]} -O ${out_file} -d ${out_dtbl} -a hmmer_alignment -s ${sfld_params.last()[1]} -o hmmer_sfld_processed
    else
        # If the file doesn't exist, create an empty file called "hmmer_sfld_processed"
        # this due to no significant hits being found in the batch
        touch "hmmer_sfld_processed"
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
