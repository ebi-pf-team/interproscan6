include { CHECK_NUCLEIC } from "$projectDir/interproscan/modules/pre_checks/main"
include { CHECK_DATA } from "$projectDir/interproscan/subworkflows/sequence_analysis/check_data"
include { CHECK_XREF_DATA } from "$projectDir/interproscan/subworkflows/xrefs/check_xref_data"

def printHelp() {
    """
    Usage example:
        nextflow run interproscan.nf --input <path to fasta file>

    Params options:
        --applications <ANALYSES>          Optional, comma separated - without spaces - list of analysis methods
                                            (i.e. member databases/applications).
                                            If this option is not set, ALL Interpro consortium member analyses will be run.
        --datadir                          Optional, path to the data dir. Default 'data' in the Interproscan 
                                            project dir.
        --disable_precalc                  Optional. Disables use of the precalculated match lookup service.
                                            All match calculations will be run locally.
        --formats <FORMATS> Optional, comma separated - without spaces - list of output formats.
                                            Accepted: tsv, json and xml
        --goterms Optional. Include GO terms in the output.
        --help                             Optional, display help information
        --input <INPUT-FILE-PATH>          [REQUIRED] Path to fasta file that should be loaded on Master startup.
        --nucleic                          Optional. Input comprises nucleic acid sequences. [Boolean]
        --outdir <OUTPUT-DIR-PATH>         Optional. Path to the output dir.
                                            Output files are automatically named after the input file, with the 
                                            suffix '.ips6.*'. Default: present working dir.
        --pathways Optional. Include pathway information in the output. [Boolean]
        --signalp_mode Optional. Set which SignalP/SignalP_EUK prediction models are used. Models may have to be installed.
                                            Accepted: 'fast', 'slow', 'slow-sequential'. Default: 'fast'.
        --version                          Print the version of InterProScan.
    """
}

workflow PRE_CHECKS {
    take:
    help_msg
    seq_input
    data_dir
    using_nucleic
    all_params
    user_applications
    output_formats
    version_msg
    ipscn_version
    signalp_mode
    signalp_gpu
    goterms
    pathways

    main:
    if ( !nextflow.version.matches('>=23.04') ) {
        println "InterProScan requires Nextflow version 23.04 or greater -- You are running version $nextflow.version"
        exit 1
    }

    if (help_msg) {
        log.info printHelp()
        exit 0
    }

    if (version_msg) {
        log.info "InterProScan version: ${ipscn_version}"
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
    def parameters_expected = [
        'input', 'applications', 'disable_precalc', 'help', 'datadir',
        'batchsize', 'url_precalc', 'check_precalc', 'matches',
        'sites', 'bin', 'members', 'translate', 'nucleic',
        'formats', 'outdir', 'xrefs', 'goterms', 'pathways', 'signalp_mode',
        'ipscn_version', 'version', 'lookup_retries', 'signalp_gpu'
    ]
    def parameter_diff = all_params - parameters_expected
    if (parameter_diff.size() != 0) {
        log.info printHelp()
        exit 22, "Input not valid: $parameter_diff"
    }

    // Check if the applications are valid
    def applications_expected = [
        'antifam', 'cdd', 'coils', 'funfam', 'gene3d', 'hamap',
        'mobidb', 'ncbifam', 'panther', 'pfam', 'phobius','pirsf', 'pirsr',
        'prints', 'prosite_patterns', 'prosite_profiles',
        'sfld', 'signalp', 'signalp_euk', 'smart', 'superfamily'
    ]

    def applications_diff = user_applications.toLowerCase().split(',') - applications_expected
    if (applications_diff.size() != 0){
        log.info printHelp()
        exit 22, "Applications not valid: $applications_diff. Valid applications are: $applications_expected"
    }

    if ("${signalp_mode}".toLowerCase() !in ['fast', 'slow', 'slow-sequential']) {
        log.error "Unrecognised SignalP mode '${signalp_mode}'.\nAccepted modes: 'fast', 'slow', 'slow-sequential'"
        exit 22
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

    applications = user_applications.toLowerCase()

    CHECK_DATA(applications, data_dir)
    missingData = CHECK_DATA.out.missingData.val
    dataDir = CHECK_DATA.out.dataDir.val

    CHECK_XREF_DATA(dataDir, goterms, pathways)
    missingXrefs = CHECK_XREF_DATA.out.missingXrefs.val
    xrefDataDir = CHECK_XREF_DATA.out.xrefDataDir.val
    // if the triple quoted text has 0 indents the indents are included
    // in the output
    if (missingData && missingXrefs) {
        log.error """
Could not find all necessary data files in '${dataDir}/' and xref files in '${dataDir}/${xrefDataDir}/'
Missing files:
${missingData}
${missingXrefs}
        """
        exit 5
    } else if (missingData && !missingXrefs) {
        log.error """
Could not find all necessary data files in '${dataDir}/'
Missing files:
${missingData}
        """
        exit 5
    } else if (missingXrefs && !missingData) {
        log.error """
Could not find all necessary XREF files in '${dataDir}/${xrefDataDir}/'
Missing files:
${missingXrefs}
        """
        exit 5
    }

    log.info "Number of sequences to analyse: ${seq_input.countFasta()}"

    emit:
    dataDir
}
