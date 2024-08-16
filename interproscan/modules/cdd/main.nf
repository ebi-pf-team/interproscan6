process CDD_RUNNER {
    label 'cdd_runner'

    input:
    tuple path(fasta), val(library), val(switches), val(postprocessing_params)

    output:
    path "rpsblast_out"
    val postprocessing_params

    script:
    """
    export LD_LIBRARY_PATH="/opt/blast/lib"
    /opt/blast/rpsblast -query ${fasta} -db ${library} -out rpsblast_out ${switches}
    """
}


process CDD_POSTPROCESS {
    /*
    Process output from offline rpsblast to annotate sequence with conserved
    domain information

    rpsbproc desc:
    This utility processes domain hit data produced by local RPS-BLAST and
    generates domain family and/or superfamily annotations on the query
    sequences. Instead of retrieving domain data from database, this program
    processes dumped datafiles to obtain required information. All data files
    are downloadable from NCBI ftp site. Read README file for details
    */
    label 'analysis_parser'

    input:
    path "rpsblast_out"
    val postprocessing_params
    /*
    [0] = switches
    [1] = data
    */

    output:
    path "rpsblast_processed"
    script:
    """
    /opt/rpsbproc/RpsbProc-x64-linux/rpsbproc --infile rpsblast_out --outfile rpsblast_processed ${postprocessing_params[0]} --data-path ${postprocessing_params[1]}
    """
}


process CDD_PARSER {
    label 'analysis_parser'

    input:
    path rpsblast_processed

    output:
    path "cdd_parsed.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/cdd/cdd_parser.py \\
        ${rpsblast_processed} \\
        cdd_parsed.json
}
