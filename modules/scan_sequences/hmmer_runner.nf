process HMMER_RUNNER {
    container 'interproscan6'

    input:
    tuple path(sequences), path(hmm), val(switches)

    output:
    path "hmmer_${hmm}_${sequences}.out"
    path "hmmer_${hmm}_${sequences}.tbl"

    script:
    """
    ${params.hmmsearch_bin} ${switches} -o hmmer_${hmm}_${sequences}.out --tblout hmmer_${hmm}_${sequences}.tbl --domtblout hmmer_${hmm}_${sequences}.dtbl ${hmm} ${sequences}
    """
}
