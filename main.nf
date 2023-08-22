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
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
// def valid_params = [
//     input         : ['*.fasta'],
//     applications  : ['antifam', 'cdd', 'coils', 'funfam', 'gene3d', 'hamap', 'mobidblite', 'ncbifam', 'panther',
//                      'pfam', 'phobius', 'pirsf', 'pirsr', 'prints', 'prositepatterns', 'prositeprofiles', 'sfld',
//                      'signalp_euk', 'signalp_gram_negative', 'signalp_gram_positive', 'smart', 'super_family', 'tmhmm'],
//     formats       : ['tsv', 'xml', 'json', 'gff3'],
//     goterms       : [true, false],
//     pathways      : [true, false]
// ]


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MATCHLOOKUP } from "$projectDir/modules/lookup/match_lookup"
include { MAIN_SCAN } from "$projectDir/modules/scan_sequences/main_scan"
include { XREFS } from "$projectDir/modules/xrefs"
include { WRITERESULTS } from "$projectDir/modules/write_results"


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

    entries_path = params.xref.entries
    applications = []
    goterms_path = false
    pathways_path = false

    if (params.goterms) {
        goterms_path = params.xref.goterms
    }
    if (params.pathways) {
        pathways_path = params.xref.pathways
    }
    if (params.applications) {
        applications = Channel.of(params.applications)
    }

    if (params.disable_precalc){
        MAIN_SCAN(fasta_file, applications)
        input_xrefs = MAIN_SCAN.out
    }
    else{
        MATCHLOOKUP(fasta_file, applications)
        input_xrefs = MATCHLOOKUP.out
    }

//     XREFS(input_xrefs, entries_path, goterms_path, pathways_path)
//     WRITERESULTS(XREFS.out, params.formats)
}
