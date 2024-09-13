process MOBIDB_RUNNER {
    label 'mobidb_runner'

    /*
    no switches needed for idrpred for now
    */
    input:
    tuple path(fasta), val(switches), val(release)

    output:
    path "mobidb_out.tsv"
    val release

    script:
    """
    idrpred ${fasta} mobidb_out.tsv
    """
}


process MOBIDB_PARSER {
    label 'analysis_parser'

    input:
    path mobidb_out
    val release

    output:
    path "mobidb_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/mobidb/parser.py \\
        ${mobidb_out} \\
        ${release} \\
        mobidb_parsed.json
    """
}
