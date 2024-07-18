process MOBIDB_RUNNER {
    container 'docker.io/library/idrpred'
    label 'mobidb_runner'

    /*
    no switches needed for idrpred for now
    */
    input:
    tuple path(fasta), val(release), val(switches)

    output:
    path "idrpred_out.tsv"
    val release

    script:
    """
    idrpred ${fasta} idrpred_out.tsv
    """
}


process MOBIDB_FILTER {
    label 'analysis_filter'

    input:
    path idrpred_out
    val release

    output:
    path "mobidb_filtered.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/mobidb/parser.py \\
        ${idrpred_out} ${release} > mobidb_filtered.json
    """
}
