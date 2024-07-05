process TMHMM_RUNNER {
    container 'docker.io/library/tmhmm'
    label 'tmhmm_runner'

    input:
    tuple path(fasta), val(release)

    output:
    path "biolib_results"
    val release

    script:
    """
    biolib run DTU/DeepTMHMM --fasta ${fasta}
    """

}


process TMHMM_PARSER {
    label 'analysis_parser'

    input:
    path biolib_results
    val release

    output:
    path "tmhmm_parsed.json"

    script:
    """
    python3 $projectDir/scripts/tmhmm/parser.py \
        ${biolib_results}/TMRs.gff3 \
        ${release}
        > tmhmm_parsed.json
    """
}
