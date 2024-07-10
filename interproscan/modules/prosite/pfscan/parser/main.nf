process PFSCAN_PARSER {
    /*
    Parse the output from the ps_scan.pl script (wrapper for pfscan)
    into the internal IPS6 JSON structure
    */
    label 'analysis_parser'

    input:
        path ps_scan_out

    output:
        path "${ps_scan_out}_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prosite/pfscan_parser.py ${ps_scan_out} ${ps_scan_out}_parsed.json
    """
}
