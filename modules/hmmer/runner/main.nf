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
    used in the search. The results are concatenated into a single file following
    post-processing by cath-resolve-hits in a downstream process.
    */
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'
    label 'hmmer_runner'

    input:
        tuple path(fasta), path(hmm), val(switches), val(release), val(postprocessing_params), path(gene3d_cath_superfamilies)
    /*
    gene3d_cath_superfamilies is a plain text file listing all Cath superfamilies
    against which gene3D generated a hit.
    */

    output:
        path "${release}_${hmm}*"
        path "${release}_${hmm}*"
        path "${hmm}_alignment"
        val postprocessing_params

    script:
    """
    awk 'NR>1' ${gene3d_out} | while read line
    do
        item=\$(echo \$line | cut -f2 -d ' ')
        new_item=\${item//./\\/}
        hmm_file_path="\${hmm}/\${new_item}.hmm"
        if [ -f \$hmm_file_path ]; then
            hmmsearch ${switches} -o ${release}_\${hmm_file_path}.out --domtblout ${release}_\${hmm_file_path}.dtbl \$hmm_file_path ${fasta}
        fi
    done
    touch ${hmm}_alignment
    """
}
