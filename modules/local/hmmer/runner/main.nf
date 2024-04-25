process HMMER_RUNNER {
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'

    input:
    tuple path(fasta), path(hmm), val(switches), val(release), val(alignment), val(postprocessing_params)
    /*
    The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
    But only generate alignmnets for these tool to reduce volume size
    */

    output:
    path "${hmm}_${release}.out"
    path "${hmm}_${release}.dtbl"
    path "${hmm}_alignment"
    val postprocessing_params

    script:
    """
    hmmsearch ${switches} -o ${hmm}_${release}.out --domtblout ${hmm}_${release}.dtbl ${alignment ? "-A ${hmm}_alignment" : ""} ${hmm} ${fasta}
    if [ ! -f ${hmm}_alignment ]; then
        touch ${hmm}_alignment
    fi
    """
}
