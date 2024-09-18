include { CHECK_SEQUENCES } from "$projectDir/interproscan/modules/pre_checks/main"
include { CHECK_DATA } from "$projectDir/interproscan/subworkflows/sequence_analysis/check_data"
include { CHECK_XREF_DATA } from "$projectDir/interproscan/subworkflows/xrefs/check_xref_data"

def printHelp() {
    """
Usage example:
    nextflow run interproscan.nf -profile <executor, container runtime> --input <path to fasta file>

Arguments:
    [Required]
    -profile                    Define the runtime profiles to use. Note the signal dash!
                                    Define an executor (built-in: 'local', 'slurm' and 'lsf') 
                                    and container runtime (built-in: 'docker', 'singularity', and 'apptainer')
    --input <INPUT-FILE-PATH>   Path to fasta file of sequences to be analysed.


    [Optional]
    --applications <ANALYSES>   Comma separated (without spaces) listing applications/member DBs to run.
                                    Default: All Interpro consortium members (except MobiDB-Lite) will run.
    --datadir <DATA-DIR>        Path to the data dir. Default 'data' dir in the Interproscan project dir.
                                    Nextflow does not tolerate spaces in paths.
    --disable_precalc           Disables use of the precalculated match lookup service. [Boolean]
                                    All match calculations will be run locally.
    --formats <FORMATS>         Comma separated (without spaces) list output file formats.
                                    Accepted: tsv, json and xml. Default: tsv,json,xml
    --goterms                   Include GO terms in the output. [Boolean]
    --help                      Display help information. [Boolean]
    --nucleic                   Input comprises nucleic acid sequences. [Boolean]
    --outdir <OUTPUT-DIR>       Path to the output dir.
                                    Output files are automatically named after the input file, with the 
                                    suffix '.ips6.*'. Default: present working dir.
                                    Nextflow does not tolerate spaces in paths.
    --pathways                  Include pathway information in the output. [Boolean]
    --signalp_mode <MODE>       SignalP/SignalP_EUK prediction models to use. Models may have to be installed.
                                    Accepted: 'fast', 'slow', 'slow-sequential'. Default: 'fast'.
    --version                   Print the version of InterProScan.[Boolean]

Please give us your feedback by sending an email to

interhelp@ebi.ac.uk

Copyright © EMBL European Bioinformatics Institute, Hinxton, Cambridge, UK. (http://www.ebi.ac.uk) The InterProScan
software itself is provided under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0.html).
Third party components (e.g. member database binaries and models) are subject to separate licensing - please see the
individual member database websites for details.

Built-in analyses:
    AntiFam: Profile-HMMs designed to identify spurious protein predictions.
    CDD: Predict protein domains and families based on well-annotated multiple sequence alignment models.
    Coils: Prediction of coiled coil regions in proteins.
    Hamap: High-quality Automated and Manual Annotation of Microbial Proteomes.
    Gene3D: Structural assignment for whole genes and genomes using the CATH domain structure database.
    FunFam: Protein function annotations for protein families and superfamilies, based upon evolutionary relationships
    NCBIfam: NCBIFams (including the original TIGRFAMs) are protein families based on hidden Markov models (HMMs).
    PANTHER: The PANTHER (Protein ANalysis THrough Evolutionary Relationships) Classification System 
            classifies genes by their functions, using published scientific experimental evidence 
            and evolutionary relationships to predict function.
    Pfam: A large collection of protein families, each represented by multiple sequence alignments and 
            hidden Markov models (HMMs).
    PIRSF: The PIRSF concept is used as a guiding principle to provide comprehensive and 
            non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect 
            their evolutionary relationships.
    PIRSR: PIRSR is a database of protein families based on hidden Markov models (HMMs) and Site Rules.
    PRINTS: A compendium of protein fingerprints - a fingerprint is a group of conserved motifs 
            used to characterise a protein family.
    ProSite_Patterns: Documentation entries describing protein domains, families and functional sites.
            PROSITE patterns are simple, descriptive motifs representing conserved sequences.
    ProSite_Profiles: Documentation entries describing protein domains, families and functional sites.
            PROSITE profiles are detailed, position-specific scoring matrices that offer a more sensitive 
            and comprehensive means of identifying and classifying protein domains and families.
    SFLD: SFLD is a database of protein families based on hidden Markov models (HMMs).
    SUPERFAMILY: SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
    SMART: SMART allows the identification and analysis of domain architectures based on hidden Markov models (HMMs).

Licensed analyses (require additional installation steps):
    DeepTMHMM: Coming Soon!
    MobiDB: Prediction of intrinsically disordered regions in proteins. 
            Runs idrpred to check for hits against a MobiDB-Lite database.
    Phobius:  A combined transmembrane topology and signal peptide predictor.
    SignalP: Signal peptide prediction using all SignalP models.
    SignalP_EUK : Signal peptide prediction using SignalP, and triggers post-processing of the SP 
            predictions by SignalP6 to prevent spurious results (only predicts type Sec/SPI).
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
    outdir

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
    // and it always checks the input FASTA file for illegal characters
    // this includes member specific and general illegal characters
    if (using_nucleic) {
        is_nucleic = true
    } else {
        is_nucleic = false
    }
    CHECK_SEQUENCES(seq_input, seq_input, is_nucleic, outdir)

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
