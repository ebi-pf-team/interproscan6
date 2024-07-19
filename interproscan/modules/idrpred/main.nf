process IDRPRED_RUNNER {
    container 'docker.io/library/idrpred'
    label 'idrpred_runner'

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


process IDRPRED_PARSER {
    label 'analysis_parser'

    input:
    path idrpred_out
    val release

    output:
    path "idrpred_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/idrpred/parser.py \\
        ${idrpred_out} ${release} > idrpred_parsed.json
    """
}
