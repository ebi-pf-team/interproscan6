process PFSCAN_PARSER {
    /*
    Parse the output from the ps_scan.pl script (wrapper for pfscan)
    into the internal IPS6 JSON structure
    */
    label 'analysis_parser'

    input:
        path ps_scan_output

    output:
        path "ips6_ps_scan.json"

    script:
    """
    $projectDir/scripts/members/prosite/pfscan_parser.py ${ps_scan_output} > ${ips6_ps_scan.json}
    """
}