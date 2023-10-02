nextflow.enable.dsl=2
import org.yaml.snakeyaml.Yaml

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HELP MESSAGE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def helpMessage() {
  log.info """
        !!!!!!!!!! OUTOFDATE !!!!!!!!!!!!
        SEE README!

        Usage:
        The typical command for running the pipeline is as follows:
        nextflow run main.nf --input ./files_test/test_all_appl.fasta

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
include { HASH_SEQUENCE } from "$projectDir/modules/hash_sequence"
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
    entries_path = params.xref.entries
    output_path = input_yaml.outfile

    Channel.fromPath( input_yaml.input )
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { fasta_channel }

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

    HASH_SEQUENCE(fasta_channel)

    lookup_to_scan = null
    matches_lookup = []
    if (!input_yaml.disable_precalc) {
        MATCHLOOKUP(HASH_SEQUENCE.out, applications)
        matches_lookup = MATCHLOOKUP.out.map { it.first() }
        lookup_to_scan = MATCHLOOKUP.out.map { it.last() }
    }

    if (input_yaml.disable_precalc || lookup_to_scan) {
        applications_channel = Channel.fromList(applications)
        if (lookup_to_scan) {
            fasta_application = lookup_to_scan.combine(applications_channel)
        }
        else {
            fasta_application = fasta_channel
            .combine(applications_channel)
        }
        MAIN_SCAN(fasta_application)
    }

    // I need to improve matches_lookup output and join it with MAIN_SCAN.out before XREFS!!
    XREFS(MAIN_SCAN.out, entries_path, goterms_path, pathways_path)

    Channel.fromList(input_yaml.formats)
    .set { formats_channel }

    XREFS.out
    .collect()
    .set { collected_outputs }

    HASH_SEQUENCE.out
    .collect()
    .set { collected_sequences }

    WRITERESULTS(collected_sequences, collected_outputs, formats_channel, output_path)
}





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