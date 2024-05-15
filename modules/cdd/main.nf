process CDD_RUNNER {
    label 'cdd_runner'

    container 'docker.io/ncbi/blast'

    input:
    tuple path(fasta), val(library), val(release), val(switches), val(postprocessing_params)

    output:
    path "rpsblast_out"
    val release
    val postprocessing_params

    script:
    """
    rpsblast -query ${fasta} -db ${library} -out rpsblast_out ${switches}
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
    container 'docker.io/library/cdd'
    label 'analysis_parser'

    input:
    path "rpsblast_out"
    val release
    val postprocessing_params
    /*
    [0] = switches
    [1] = data
    */

    output:
    path "rpsblast_processed"
    val release
    script:
    """
    /opt/cdd/RpsbProc-x64-linux/rpsbproc --infile rpsblast_out --outfile rpsblast_processed ${postprocessing_params[0]} --data-path ${postprocessing_params[1]}
    """
}


process CDD_PARSER {
    label 'analysis_parser'

    input:
    path rpsblast_processed
    val release

    output:
    path "cdd_parsed.json"

    script:
    """
    python3 $projectDir/scripts/members/cdd/cdd_parser.py \\
        ${rpsblast_processed} \\
        ${release} > cdd_parsed.json
    """
}
