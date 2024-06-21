process HMMER_RUNNER {
    label 'hmmer_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(release), val(postprocessing_params)

    output:
        path "${release}._.${member}._.out"
        val postprocessing_params
    script:
    """
    /opt/hmmer/bin/hmmsearch ${switches} -o ${release}._.${member}._.out ${hmm} ${fasta}
    """
}

process HMMER_RUNNER_WITH_ALIGNMENTS {
    label 'hmmer_runner'
    /*
    The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
    */
    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(release), val(postprocessing_params)

    output:
        path "${release}._.${member}._.*"
        val postprocessing_params
        path "${member}_alignment"

    script:
    """
    /opt/hmmer/bin/hmmsearch ${switches} \
    -o ${release}._.${member}._.out \
    --domtblout ${release}._.${member}._.dtbl \
    -A ${member}_alignment" ${hmm} ${fasta}
    """
}


process HMMER_RUNNER_TBL_OUTPUT {
    label 'hmmer_runner'
    /*
    The post processing of HAMAP requires the tbl file
    */
    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(release), val(postprocessing_params)

    output:
        path "${release}._.${member}._.out"
        val postprocessing_params
        path "${release}._.${member}._.table.tbl"

    script:
    """
    /opt/hmmer/bin/hmmsearch ${switches} -o ${release}._.${member}._.out --tblout ${release}._.${member}._.table.tbl ${hmm} ${fasta}
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
    label 'hmmer_runner'

    when:
        "${applications}".contains('funfam')

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(release), val(build_alignment), val(build_table), val(postprocessing_params), val(cath_superfamily)
        val applications

    /*
    post-processing params:
    4. FunFam HMM dir
    5. FunFam HMMsearch switches
    6. FunFam release number
    */

    output:
        path "${postprocessing_params[6]}._.funfam._.${cath_superfamily}.out"
        path "${postprocessing_params[6]}._.funfam._.${cath_superfamily}.dtbl"
        val postprocessing_params

    script:
    """
    /opt/hmmer/bin/hmmsearch \\
        ${postprocessing_params[5]} \\
        -o ${postprocessing_params[6]}._.funfam._.${cath_superfamily}.out \\
        --domtblout ${postprocessing_params[6]}._.funfam._.${cath_superfamily}.dtbl \\
        "${postprocessing_params[4]}${cath_superfamily.replace('.', '/')}.hmm" \\
        ${fasta}
    """

}
