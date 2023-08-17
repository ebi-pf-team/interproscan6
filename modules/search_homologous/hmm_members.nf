process MEMBERS_HMM {
    input:
    val fasta
    val applications

    output:
    path 'parsed_hmm'

    script:
    """
    python $projectDir/scripts/homologs_search/hmm_members.py -fasta ${fasta} -appl ${applications} > parsed_hmm
    """
}