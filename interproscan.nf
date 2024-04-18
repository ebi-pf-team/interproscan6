nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"
include { GET_ORFS } from "$projectDir/modules/local/get_orfs/main"
include { SEQUENCE_PRECALC } from "$projectDir/subworkflows/sequence_precalc/main"
include { SEQUENCE_ANALYSIS } from "$projectDir/subworkflows/sequence_analysis/main"
include { AGGREGATE_RESULTS } from "$projectDir/modules/local/write_output/aggregate_results/main"


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
        --seqtype <SEQUENCE-TYPE>          Optional. Optional, the type of the input sequences (dna/rna (n) or protein (p)). The default sequence type is protein.
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
            The typical command for running the pipeline is:
                nextflow run interproscan.nf --input <path to fasta file>
            For more information, please use the --help flag.
            """
            exit 1
}

// Check if the input parameters are valid
def parameters_expected = ['input', 'applications', 'disable_precalc', 'help', 'batchsize', 'url_precalc', 'check_precalc', 'matches', 'sites', 'bin', 'members', 'tsv_pro', 'translate', 'nucleic', 'orfs']
def parameter_diff = params.keySet() - parameters_expected
if (parameter_diff.size() != 0){
    log.info printHelp()
    exit 1, "Input not valid: $parameter_diff"
}

// Check if the applications are valid
def applications_expected = ['antifam', 'cdd', 'coils', 'funfam', 'gene3d', 'hamap', 'mobidblite', 'ncbifam', 'panther', 'pfam', 'phobius', 'pirsf', 'pirsr', 'prints', 'prositepatterns', 'prositeprofiles', 'sfld', 'signalp', 'smart', 'superfamily', 'tmhmm']
def applications_diff = params.applications.toLowerCase().split(',') - applications_expected
if (applications_diff.size() != 0){
    log.info printHelp()
    exit 1, "Applications not valid: $applications_diff. Valid applications are: $applications_expected"
}

// Check if the input file is a fasta file and if it contains sequences
if (!params.input.toLowerCase().find(/.fasta$|.faa$|.fna$/)) {
    log.error "The input file is not a FASTA file (it does not end in .fasta, .faa or .fna)"
    exit 1
}

def seq_count = file(params.input).countFasta()
    if (seq_count == 0) {
        log.info "No sequence found in the input file"
        exit 1
    }

log.info "Number of sequences to analyse: ${seq_count}"

workflow {
    applications = params.applications.toLowerCase()

    Channel.fromPath( params.input , checkIfExists: true)
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { ch_fasta }

    if (params.nucleic) {
        if (params.translate.strand.toLowerCase() !in ['both','plus','minus']) {
            log.info "Strand option '${params.translate.strand.toLowerCase()}' in nextflow.config not recognised. Accepted: 'both', 'plus', 'minus'"
            exit 1
        }
        GET_ORFS(ch_fasta, params.translate.strand, params.translate.methionine, params.translate.min_len, params.translate.genetic_code)
        GET_ORFS.out.splitFasta( by: params.batchsize, file: true )
        .set { orfs_fasta }
        PARSE_SEQUENCE(orfs_fasta)
    }
    else {
        PARSE_SEQUENCE(ch_fasta)
    }

    sequences_to_analyse = null
    parsed_matches = Channel.empty()
    if (!params.disable_precalc) {
        log.info "Using precalculated match lookup service"
        SEQUENCE_PRECALC(PARSE_SEQUENCE.out, applications)
        parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
        sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
    }

    analysis_result = Channel.empty()
    if (params.disable_precalc || sequences_to_analyse) {
        log.info "Running sequence analysis"
        if (sequences_to_analyse) {
            fasta_to_runner = sequences_to_analyse
        }
        else {
            fasta_to_runner = ch_fasta
        }
        analysis_result = SEQUENCE_ANALYSIS(fasta_to_runner, applications)
    }

    all_results = parsed_matches.collect().concat(analysis_result.collect())

    AGGREGATE_RESULTS(all_results.collect())
    AGGREGATE_RESULTS.out.view()
}
