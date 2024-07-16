process MOBIDB_RUNNER {
    label 'mobidb_runner'

    input:
    tuple path(fasta), val(library), val(release), val(switches), val(postprocessing_params)

    output:
    path "idrpred_out"
    val release
    val postprocessing_params

    script:
    """
    /opt/mobidb/idrpred ${switches} --infile ${fasta} --outfile idrpred_out
    """
}


process MOBIDB_FILTER {
    label 'analysis_filter'

    input:
    path idrpred_out

    output:
    path "mobidb_filtered.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/mobidb/mobidb_filter.py \\
        ${idrpred_out} > mobidb_filtered.json
    """
}
