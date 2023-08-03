nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
  log.info """
        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --input ./files_test/test_all_appl.fasta -params-file input_opt.yaml

        Mandatory arguments:
        --input                     Directory to input fasta file

        Optional arguments:
        --help                      This usage statement.
        --input-opt                 YAML file with analyses and output options
        --goterms                   Switch on lookup of corresponding Gene Ontology.
        --pathways                  Switch on lookup of corresponding Pathway.
        """
}

if (params.help) {
    helpMessage()
    exit 0
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MATCHLOOKUP } from "$projectDir/modules/lookup/match_lookup"
include { XREFS } from "$projectDir/modules/xrefs"
include { WRITERESULTS } from "$projectDir/modules/lookup/write_results"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow {
    ch_input = file(params.input, checkIfExists: true)
    if (ch_input.isEmpty()) {
        exit 1, "File provided with --input is empty: ${ch_input.getName()}!"
    }
    Channel.fromPath( ch_input )
    .splitFasta( by: params.batchsize, file: true )
    .set { fasta_file }

    entries_path = params.data.entries

    if (params.goterms) {
        goterms_path = params.data.goterms
    }
    if (params.pathways) {
        pathways_path = params.data.pathways
    }

    MATCHLOOKUP(fasta_file, params.applications)
    XREFS(MATCHLOOKUP.out, entries_path, goterms_path, pathways_path)
//     WRITERESULTS(XREFS.out, params.formats)
}
