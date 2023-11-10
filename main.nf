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
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"
include { UNION_RESULTS } from "$projectDir/modules/local/union_results/main"
include { XREFS } from "$projectDir/modules/local/xrefs/main"
include { WRITE_RESULTS } from "$projectDir/modules/local/write_results/main"
include { SEQUENCE_PRECALC } from "$projectDir/subworkflows/sequence_precalc/main"
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

    Channel.fromPath( input_yaml.input )
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { fasta_channel }

    PARSE_SEQUENCE(fasta_channel)

    sequences_to_analyse = null
    parsed_matches = null
    if (!input_yaml.disable_precalc) {
        SEQUENCE_PRECALC(PARSE_SEQUENCE.out, applications)
        parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
        sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
    }

    analysis_result = null
    if (input_yaml.disable_precalc || sequences_to_analyse) {
        applications_channel = Channel.fromList(applications)
        if (sequences_to_analyse) {
            fasta_application = sequences_to_analyse.combine(applications_channel)
        }
        else {
            fasta_application = fasta_channel
            .combine(applications_channel)
        }
        analysis_result = SEQUENCE_ANALYSIS(fasta_application, check_tsv_pro)
    }

    if (parsed_matches) {
        parsed_matches
        .collect()
        .set { all_parsed_lookup }
    }
    else {
        all_parsed_lookup = []
    }

    if (analysis_result) {
        analysis_result
        .collect()
        .set { all_parsed_analysis }
    }
    else {
        all_parsed_analysis = []
    }

    UNION_RESULTS(all_parsed_lookup, all_parsed_analysis)
    XREFS(UNION_RESULTS.out, entries_path, goterms_path, pathways_path)

    PARSE_SEQUENCE.out
    .collect()
    .set { all_sequences_parsed }

    Channel.fromList(input_yaml.formats)
    .set { formats_channel }

    WRITE_RESULTS(all_sequences_parsed, XREFS.out, formats_channel, output_path)
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
