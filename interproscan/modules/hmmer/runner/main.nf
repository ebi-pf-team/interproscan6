process HMMER_RUNNER {
    label 'hmmer_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)

    output:
        path "${member}_out"
        val postprocessing_params
        val member

    script:
    """
    /opt/hmmer3/bin/hmmsearch ${switches} -o ${member}_out ${hmm} ${fasta}
    """
}


process HMMER_RUNNER_WITH_ALIGNMENTS {
    label 'hmmer_runner'
    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)

    output:
        path "${member}_out"
        val postprocessing_params
        val member
        path "${member}_alignment"
        path "${member}_dtbl"

    script:
    """
    /opt/hmmer3/bin/hmmsearch ${switches} -o ${member}_out --domtblout ${member}_dtbl -A ${member}_alignment ${hmm} ${fasta}
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

    This process will only run in at least one Cath superfamily is specified.
    */
    label 'hmmer_runner'

    when:
        cath_superfamilies.size() > 0 && "${applications}".contains('funfam')

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)
        path cath_superfamilies
        val applications
    /*
    post-processing params:
    4. FunFam HMM dir
    5. FunFam HMMsearch switches
    */

    output:
        path "funfam_*.out"
        val postprocessing_params
        val member

    script:
    """
    while IFS= read -r cath_superfamily
    do
        hmm_file="\${cath_superfamily//./\\/}.hmm"
        /opt/hmmer3/bin/hmmsearch \\
            ${postprocessing_params[5]} \\
            -o funfam_\${cath_superfamily}.out \\
            "${postprocessing_params[4]}/\$hmm_file" \\
            ${fasta}
    done < ${cath_superfamilies}
    """
}


process HAMAP_HMMER_RUNNER {
    label 'hmmer_runner'
    /*
    The post processing of HAMAP requires the tbl file
    */
    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)

    output:
        path "${member}_out"
        val postprocessing_params
        val member
        path fasta
        path "${member}_table.tbl"

    script:
    """
    /opt/hmmer3/bin/hmmsearch \\
        ${switches} -o ${member}_out \\
        --tblout ${member}_table.tbl \\
        ${hmm} \\
        ${fasta}
    """
}


process PIRSF_HMMER_RUNNER {
    label 'hmmer_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)
    /*
    No -Z number was provided in i5, so migrating from hmmscan to hmmsearch
    results in a chnage in the E-values, so we have to keep with hmmscan
    for now.
    */

    output:
        path "${member}_out"
        val postprocessing_params
        val member
        path "${member}_dtbl"

    script:
    """
    /opt/hmmer3/bin/hmmpress ${hmm}
    /opt/hmmer3/bin/hmmscan \\
        ${switches} -o \\
        ${member}_out \\
        --domtblout ${member}_dtbl \\
        ${hmm} \\
        ${fasta}
    """
}


process SMART_HMMER2_RUNNER {
    label 'hmmer_2_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches)

    output:
        path "${member}_out"
        path "${fasta}"  // Used for filtering kinase hits in SMART_FILTER_MATCHES
        val member

    script:
    """
    /opt/hmmer2/bin/hmmpfam ${switches} ${hmm} ${fasta} > ${member}_out
    """
}


process HMMER_SCAN_RUNNER {
    label 'hmmer_runner'

    input:
        tuple path(fasta), val(member), path(hmm), val(switches), val(postprocessing_params)
    /*
    Superfamily uses a .pl script that create assignments from the output of HMMER3 hmmscan
    */

    output:
        path "${member}_out"
        val postprocessing_params
        val member
        path fasta
        path hmm

    script:
    """
    /opt/hmmer3/bin/hmmpress ${hmm}
    /opt/hmmer3/bin/hmmscan ${switches} -o ${member}_out ${hmm} ${fasta}
    """
}
