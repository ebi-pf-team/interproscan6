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
include { HASH_SEQUENCE } from "$projectDir/modules/hash_sequence/main"
include { MATCH_LOOKUP } from "$projectDir/modules/match_lookup/main"
include { XREFS } from "$projectDir/modules/xrefs/main"
include { WRITE_RESULTS } from "$projectDir/modules/write_results/main"
include { SEQUENCE_ANALYSIS } from "$projectDir/subworkflows/sequence_analysis/main"


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

    check_tsv_pro = input_yaml.formats.contains("TSV-PRO")

    HASH_SEQUENCE(fasta_channel)

    lookup_to_scan = null
    matches_lookup = []
    if (!input_yaml.disable_precalc) {
        MATCH_LOOKUP(HASH_SEQUENCE.out, applications)
        matches_lookup = MATCH_LOOKUP.out.map { it.first() }
        lookup_to_scan = MATCH_LOOKUP.out.map { it.last() }
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
        SEQUENCE_ANALYSIS(fasta_application, check_tsv_pro)
    }

    // I need to improve matches_lookup output and join it with MAIN_SCAN.out before XREFS!!
    XREFS(SEQUENCE_ANALYSIS.out, entries_path, goterms_path, pathways_path)

    Channel.fromList(input_yaml.formats)
    .set { formats_channel }

    XREFS.out
    .collect()
    .set { collected_outputs }

    HASH_SEQUENCE.out
    .collect()
    .set { collected_sequences }

    WRITE_RESULTS(collected_sequences, collected_outputs, formats_channel, output_path)
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
