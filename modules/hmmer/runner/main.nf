process HMMER_RUNNER {
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
    label 'hmmer_runner'

    input:
        tuple path(fasta), path(hmm), val(switches), val(release), val(alignment), val(postprocessing_params)
    /*
    The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
    But only generate alignmnets for these tool to reduce volume size
    */

    output:
        path "${release}_${hmm}.out"
        path "${release}_${hmm}.dtbl"
        path "${hmm}_alignment"
        val postprocessing_params

    script:
    """
    hmmsearch ${switches} -o ${release}_${hmm}.out --domtblout ${release}_${hmm}.dtbl ${alignment ? "-A ${hmm}_alignment" : ""} ${hmm} ${fasta}
    if [ ! -f ${hmm}_alignment ]; then
        touch ${hmm}_alignment
    fi
    """
}


process FUNFAM_HMMER_RUNNER {
    /*
    FunFam requires its own runner in order to only use those FunFam families
    that are associated to Cath-Gene3D superfamilies where hits were found
    (Gene3D must be run before FunFam). Otherwise running HMMER3 for all
    FunFam hmm profiles would take an extremely long time.

    There will be one hmmer.out and one hmmer.dtbl file per FunFam hmm profile
    used in the search.
    */
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
    label 'hmmer_runner'

    when:
        "${applications}".contains('funfam')

    input:
        tuple path(fasta), path(hmm), val(switches), val(release), val(alignment), val(postprocessing_params), val(cath_superfamily)
        val applications

    /*
    post-processing params:
    4. FunFam HMM dir
    5. FunFam HMMsearch switches
    6. FunFam release number
    */

    output:
        path "${release}_funfam_${cath_superfamily}.out"
        path "${release}_funfam_${cath_superfamily}.dtbl"
        path "${hmm}_alignment"
        val postprocessing_params

    script:
    """
    hmmsearch \\
        ${postprocessing_params[5]} \\
        -o ${postprocessing_params[6]}_funfam_${cath_superfamily}.out \\
        --domtblout ${postprocessing_params[6]}_funfam_${cath_superfamily}.dtbl \\
        "${postprocessing_params[4]}${cath_superfamily.replace('.', '/')}.hmm" \\
        ${fasta}

    touch ${hmm}_alignment
    """

}
