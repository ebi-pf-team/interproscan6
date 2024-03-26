nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"
include { SEQUENCE_PRECALC } from "$projectDir/subworkflows/sequence_precalc/main"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def printHelp() {
    """
    Usage example:
        nextflow run interproscan.nf --input <path to fasta file>

    Params options:
        --applications <ANALYSES>          Optional, comma separated - without spaces - list of analysis methods (i.e. member databases/applications).
                                            If this option is not set, ALL analyses will be run.
        --disable-precalc                  Optional. Disables use of the precalculated match lookup service.
                                            All match calculations will be run locally.
        --help                             Optional, display help information
        --input <INPUT-FILE-PATH>          [REQUIRED] Path to fasta file that should be loaded on Master startup.
    """
}

if (params.help) {
    log.info printHelp()
    System.exit(0)
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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

    sequences_to_analyse = null
    parsed_matches = null

    if (!params.disable_precalc) {
        log.info "Using precalculated match lookup service"
        SEQUENCE_PRECALC(PARSE_SEQUENCE.out, params.applications)
        parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
        sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
    }

    //  Just temporary to see in which folders are the partial results
    log.info "Applications: ${params.applications}"
    SEQUENCE_PRECALC.out.parsed_matches.view()
    SEQUENCE_PRECALC.out.sequences_to_analyse.view()
}
