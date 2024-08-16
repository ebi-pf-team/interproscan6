process MOBIDB_RUNNER {
    container 'docker.io/library/idrpred'
    label 'mobidb_runner'

    /*
    no switches needed for idrpred for now
    */
    input:
    tuple path(fasta), val(switches)

    output:
    path "mobidb_out.tsv"

    script:
    """
    idrpred ${fasta} mobidb_out.tsv
    """
}


process MOBIDB_PARSER {
    label 'analysis_parser'

    input:
    path mobidb_out

    output:
    path "mobidb_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/mobidb/parser.py \\
        ${mobidb_out} \\
        mobidb_parsed.json
    """
}
