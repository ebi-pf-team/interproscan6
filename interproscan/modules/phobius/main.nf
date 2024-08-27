process PHOBIUS_RUNNER {
    container 'docker.io/library/phobius'
    label 'phobius_runner'

    input: tuple path(fasta), val(release)

    output:
    path "*._.phobius_out.txt"

    script:
    """
    phobius ${fasta} > ${release}._.phobius_out.txt
    """
}


process PHOBIUS_PARSER {
    label 'analysis_parser'

    input:
    path phobius_out

    output:
    path "phobius_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/phobius/parser.py \\
        ${phobius_out} \\
        phobius_parsed.json
    """
}