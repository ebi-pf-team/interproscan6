process RUN_RPSBLAST {
    label 'cdd_runner'

    input:
    tuple val(meta), path(fasta)
    val library

    output:
    tuple val(meta), path("rpsblast.out")

    script:
    """
    export LD_LIBRARY_PATH="/opt/blast/lib"
    /opt/blast/rpsblast \
        -query ${fasta} \
        -db ${library} \
        -out rpsblast.out \
        -evalue 0.01 -seg no -outfmt 11
    """
}


process RUN_RPSPROC {
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
    tuple val(meta), val(rpsblast_out)
    path library

    output:
    tuple val(meta), path("rpsbproc.out")

    script:
    """
    /opt/rpsbproc/RpsbProc-x64-linux/rpsbproc \
        --infile ${rpsblast_out} \
        --outfile rpsbproc.out \
        --data-path ${library} \
        -m std
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
    """
}
