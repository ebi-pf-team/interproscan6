def printHelp() {
    """
    Usage example:
        nextflow run interproscan.nf --input <path to fasta file>

    Params options:
        --applications <ANALYSES>          Optional, comma separated - without spaces - list of analysis methods (i.e. member databases/applications).
                                            If this option is not set, ALL analyses will be run.
        --disable-precalc                  Optional. Disables use of the precalculated match lookup service.
                                            All match calculations will be run locally.
        --formats <FORMATS> Optional, comma separated - without spaces - list of output formats.
        --help                             Optional, display help information
        --input <INPUT-FILE-PATH>          [REQUIRED] Path to fasta file that should be loaded on Master startup.
        --nucleic                          Optional. Input comprises nucleic acid sequences.
        --output <OUTPUT-FILE-PATH>        Optional. Path to the output file.
                                            If this option is not set, the output will be write on results/ folder.
    """
}


workflow PRE_CHECKS {
    take:
    params

    main:
    if ( !nextflow.version.matches('23.10+') ) {
        println "InterProScan requires Nextflow version 23.10 or greater -- You are running version $nextflow.version"
        exit 1
    }

    if (params.help) {
        log.info printHelp()
        exit 1
    }

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
    def parameters_expected = ['input', 'applications', 'disable_precalc', 'help', 'batchsize', 'url_precalc', 'check_precalc', 'matches', 'sites', 'bin', 'members', 'tsv_pro', 'translate', 'nucleic', 'formats', 'output']
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
    def seq_count = file(params.input).countFasta()
        if (seq_count == 0) {
            log.info "No sequence found in the input file"
            exit 1
        }

    log.info "Number of sequences to analyse: ${seq_count}"
}
