process HMMER_RUNNER {
    container 'docker.io/biocontainers/hmmer:v3.2.1dfsg-1-deb_cv1'

    input:
    tuple path(fasta), path(hmm), val(switches), val(alignment), val(sfld_postprocessing_params)
    /*
    The post processing of SFLD, FunFam and Gene3D HMMER hits requires the alignment file
    But only generate alignmnets for these tool to reduce volume size
    */

    output:
    path "hmmer_${hmm}.out"
    path "hmmer_${hmm}.dtbl"

    script:
    """
    hmmsearch ${switches} -o hmmer_${hmm}.out --domtblout hmmer_${hmm}.dtbl ${alignment ? "-A hmmer_alignment" : ""} ${hmm} ${fasta}
    """
}
