process PFSCAN_RUNNER {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the
    pftools developers. It automates running pfscan for all provided patterns and
    includes post-processing of the hits.
    */
    label 'prosite_pfscan_runner'

    input:
        tuple path(fasta), path(data), path(evaluator), val(switches)
    /*
    ps_scan_params params:
    0. patterns data file dir
    1. evaluator data file
    2. ps_scan switches
    */

    output:
        path "ps_scan.out"

    script:
    """
        perl /opt/pftools/ps_scan.pl \
        ${fasta} \
        -d ${data} \
        --pfscan /opt/pftools/pfscanV3 \
        -b ${evaluator} \
        ${switches} > ps_scan.out
    """
}
