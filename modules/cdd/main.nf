process CDD_RUNNER {
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

    input:
    path "rpsblast_out"
    val release
    val postprocessing_params
    /*
    [0] = bin file
    [1] = switches
    [2] = data
    [3] = signature list
    */

    output:
    path "rpsblast_processed"
    val release
    val postprocessing_params

    script:
    """
    ${postprocessing_params[0]} --infile rpsblast_out --outfile rpsblast_processed ${postprocessing_params[1]} --data-path ${postprocessing_params[2]}
    """
}


process CDD_ADD_SIGNATURES {
    /*
    Match hits to the signature library
    */

    input:
    path "rpsblast_processed"
    val release
    val postprocessing_params // [3] = signature list

    output:
    path "rpsblast_parsed_processed"
    val release

    script:
    """
    ${postprocessing_params[0]} --infile rpsblast_out --outfile rpsblast_processed ${postprocessing_params[1]} --data-path ${postprocessing_params[2]}
    """
}



process CDD_PARSER {
    input:
    path "rpsblast_processed"
    val release
    val something

    script:
    """
    asdfghjkl
    """
}
