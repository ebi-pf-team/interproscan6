process PHOBIUS_RUNNER {
    label 'phobius_runner'

    input: path fasta

    output:
    path "phobius_out.txt"

    script:
    """
    phobius ${fasta} > phobius_out.txt
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
