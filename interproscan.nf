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

if (!params.input) {
    log.info """
            Please provide an input file.
            The typical command for running the pipeline is as follows:
            nextflow run interproscan.nf --input <path to fasta file>
            """
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
