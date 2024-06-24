process PS_SCAN_PARSER {
    label 'analysis_parser'

    input:
        path ps_scan_output

    output:
        path "ips6_ps_scan.json"

    script:
    """
    $projectDir/scripts/members/prosite/ps_scan_parser.py ${ps_scan_output} > ${ips6_ps_scan.json}
    """
}