nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def printHelp() {
    """
    Usage example:
        nextflow run interproscan.nf --input <path to fasta file>

    Params options:
        --applications <ANALYSES>          Optional, comma separated list of analyses. If this option is not set,
                                            ALL analyses will be run.
        --help                             Optional, display help information
        --input <INPUT-FILE-PATH>          Path to fasta file that should be loaded on Master startup.
    """
}

if (params.help) {
    log.info printHelp()
    System.exit(0)
}

if (!params.input) {
    log.info """
            Please provide an input file.

            The typical command for running the pipeline is as follows:
                nextflow run interproscan.nf --input <path to fasta file>

            For more information, please use the --help flag.
            """
            exit 1
}

workflow {
    Channel.fromPath( params.input , checkIfExists: true)
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { fasta_channel }

    PARSE_SEQUENCE(fasta_channel)

//  Just temporary to see in which folders are the partial results
    PARSE_SEQUENCE.out.view()
}
