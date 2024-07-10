process PFSEARCH_PARSER {
    label 'analysis_parser'

    input:
        path pfsearch_out
        path blacklist_file

    output:
        path "${pfsearch_out}-filtered.json"

    script:
    """
    python3 $projectDir/scripts/members/prosite/pfsearch_parser.py \
        ${pfsearch_out} \
        ${pfsearch_out}-filtered.json \
        ${blacklist_file}
    """
}