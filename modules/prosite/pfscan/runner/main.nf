process PS_SCAN_RUNNER {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the 
    pftools developers. It automates running pfscan for all provided patterns and 
    includes post-processing of the hits.
    */
    container 'docker.io/sibswiss/pftools'
    label 'prosite_ps_scan_runner'

    input:
        path fasta
        path data
        path evaluator
        val release
        val switches
        val post_processing  // PROSITE profiles only
    /*
    ps_scan_params params:
    0. patterns data file dir
    1. evaluator data file
    2. ps_scan switches
    */

    output:
        path "ps_scan_${release}.out"

    script:
    """
    ps_scan.pl \
        ${fasta} \
        -d ${data} \
        --pfscan pfscanV3 \
        -b ${evaluator} \
        ${switches} > ps_scan_${release}.out
    """
}