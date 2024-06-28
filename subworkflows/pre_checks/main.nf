include { CHECK_NUCLEIC } from "$projectDir/modules/pre_checks/main"

def printHelp() {
    """
    Usage example:
        nextflow run interproscan.nf --input <path to fasta file>

    Params options:
        --applications <ANALYSES>          Optional, comma separated - without spaces - list of analysis methods (i.e. member databases/applications).
                                            If this option is not set, ALL analyses will be run.
        --disable_precalc                  Optional. Disables use of the precalculated match lookup service.
                                            All match calculations will be run locally.
        --formats <FORMATS> Optional, comma separated - without spaces - list of output formats.
        --goterms Optional. Include GO terms in the output.
        --help                             Optional, display help information
        --input <INPUT-FILE-PATH>          [REQUIRED] Path to fasta file that should be loaded on Master startup.
        --nucleic                          Optional. Input comprises nucleic acid sequences.
        --output <OUTPUT-FILE-PATH>        Optional. Path to the output file.
                                            If this option is not set, the output will be write on results/ folder.
        --pathways Optional. Include pathway information in the output.
    """
}


workflow PRE_CHECKS {
    take:
    help_msg
    seq_input
    using_nucleic
    all_params
    user_applications
    output_formats

    main:
    if ( !nextflow.version.matches('>=23.10') ) {
        println "InterProScan requires Nextflow version 23.10 or greater -- You are running version $nextflow.version"
        exit 1
    }

    if (help_msg) {
        log.info printHelp()
        exit 0
    }

    if (!seq_input) {
        log.error """
                Please provide an input file.
                The typical command for running the pipeline is:
                    nextflow run interproscan.nf --input <path to fasta file>
                For more information, please use the --help flag.
                """
                exit 5
    }

    // is user specifies the input is nucleic acid seqs
    // check the input only contains nucleic acid seqs
    if (using_nucleic) {
        try {
            CHECK_NUCLEIC(seq_input)
        } catch (all) {
            println """Error in input sequences"""
            log.error """
            The '--nucleic' flag was used, but the input FASTA file
            appears to contain at least one sequence that contains a
            non-nucleic acid residue ('A','G','C','T','*','-', case insensitive).
            Please check your input is correct.
            """
            exit 1
        }
    }

    // Check if the input parameters are valid
    def parameters_expected = ['input', 'applications', 'disable_precalc', 'help', 'batchsize', 'url_precalc', 'check_precalc', 'matches', 'sites', 'bin', 'members', 'tsv_pro', 'translate', 'nucleic', 'formats', 'output', 'xrefs', 'goterms', 'pathways', 'ipsc_version']
    def parameter_diff = all_params - parameters_expected
    if (parameter_diff.size() != 0){
        log.info printHelp()
        exit 22, "Input not valid: $parameter_diff"
    }

    // Check if the applications are valid
    def applications_expected = ['antifam', 'cdd', 'funfam', 'gene3d', 'hamap', 'ncbifam', 'panther', 'pfam', 'sfld', 'signalp']
    def applications_diff = user_applications.toLowerCase().split(',') - applications_expected
    if (applications_diff.size() != 0){
        log.info printHelp()
        exit 22, "Applications not valid: $applications_diff. Valid applications are: $applications_expected"
    }

    // Check if the formats are valid
    def formats_expected = ['json', 'tsv', 'tsv-pro', 'xml', 'gff3']
    def formats_diff = output_formats.toLowerCase().split(',') - formats_expected
    if (formats_diff.size() != 0){
        log.info printHelp()
        exit 22, "Format not valid: $formats_diff. Valid formats are: $formats_expected"
    }

    // Check if the input file is a fasta file and if it contains sequences
    if (seq_input.countFasta() == 0) {
        log.error "No sequence found in the input file"
        exit 5
    }

    log.info "Number of sequences to analyse: ${seq_input.countFasta()}"
}
