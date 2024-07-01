process HMMER_RUNNER {
    label 'hmmer_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(release), val(build_alignment), val(build_table), val(postprocessing_params)
    /*
    The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
    But only generate alignmnets for these tool to reduce volume size.
    Likewise, for the HMMER table file ).tbl)
    */

    output:
        path "${release}._.${member}._.out"
        path "${release}._.${member}._.dtbl"
        val postprocessing_params
        path "${member}_alignment"
        path "${release}._.${member}._.table.tbl"
        path "${fasta}"

    script:
    """
    /opt/hmmer/bin/hmmsearch ${switches} -o ${release}._.${member}._.out --domtblout ${release}._.${member}._.dtbl ${build_alignment ? "-A ${member}_alignment" : ""} ${build_table ? "--tblout ${release}._.${member}._.table.tbl" : ""} ${hmm} ${fasta}


    if [ ! -f ${member}_alignment ]; then
        touch ${member}_alignment
    fi
    if [ ! -f ${release}._.${member}._.table.tbl ]; then
        touch ${release}._.${member}._.table.tbl
    fi
    """
}


process funfam._.HMMER_RUNNER {
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
