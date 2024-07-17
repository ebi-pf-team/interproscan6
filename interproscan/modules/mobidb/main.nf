process MOBIDB_RUNNER {
    container 'docker.io/library/idrpred'
    label 'mobidb_runner'

    /*
    no switches needed for idrpred for now
    */
    input:
    tuple path(fasta), val(library), val(release), val(switches)

    output:
    path "idrpred_out"
    val library
    val release

    script:
    """
    idrpred ${fasta} idrpred_out
    """
}


process MOBIDB_FILTER {
    label 'analysis_filter'

    input:
    path idrpred_out
    val library
    val release

    output:
    path "mobidb_filtered.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/mobidb/mobidb_filter.py \\
        ${idrpred_out} ${library} ${release} > mobidb_filtered.json
    """
}
