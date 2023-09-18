nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

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
def f = file('input_opt.yaml')
def input_yaml = new Yaml().load(f)
def all_appl = ['AntiFam', 'CDD', 'Coils', 'FunFam', 'Gene3d', 'HAMAP', 'MobiDBLite', 'NCBIfam', 'Panther', 'Pfam',
                'Phobius', 'PIRSF', 'PIRSR', 'PRINTS', 'PrositePatterns', 'PrositeProfiles', 'SFLD', 'SignalP_EUK',
                'SignalP_GRAM_NEGATIVE', 'SignalP_GRAM_POSITIVE', 'SMART', 'SuperFamily', 'TMHMM']

workflow {
    Channel.fromPath( input_yaml.input )
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { sequences_channel }

    entries_path = params.xref.entries
    output_path = input_yaml.outfile

// SERIAL_GROUP = "PROTEIN"
//  Default for protein sequences are TSV, XML and GFF3, for nucleotide sequences GFF3 and XML.
//     if (input_yaml.formats) {
//         output_format = input_yaml.formats
//     }
//     else {
//         if SERIAL_GROUP == "PROTEIN" {
//             output_format = ["TSV", "XML", "GFF3"]
//         } else {
//             output_format = ["XML", "GFF3"]
//         }
//     }

    goterms_path = ""
    pathways_path = ""
    if (input_yaml.goterms) {
        goterms_path = params.xref.goterms
    }
    if (input_yaml.pathways) {
        pathways_path = params.xref.pathways
    }

    if (input_yaml.applications) {
        applications = input_yaml.applications
    }
    else {
        applications = all_appl
    }

    if (input_yaml.disable_precalc){
        applications_channel = Channel.fromList(applications)
        sequences_application = sequences_channel.combine(applications_channel)
        MAIN_SCAN(sequences_application)
        input_xrefs = MAIN_SCAN.out
    }
    else{
        MATCHLOOKUP(sequences_channel, applications)
        input_xrefs = MATCHLOOKUP.out
    }

    XREFS(input_xrefs, entries_path, goterms_path, pathways_path)

    XREFS.out
    .collect()
    .set { collected_outputs }

    Channel.fromList(input_yaml.formats)
    .set { formats_channel }
    WRITERESULTS(collected_outputs, formats_channel, output_path)
}
