nextflow.enable.dsl=2

include { PARSE_SEQUENCE } from "$projectDir/interproscan/modules/parse_sequence/main"
include { GET_ORFS } from "$projectDir/interproscan/modules/get_orfs/main"
include { REPRESENTATIVE_DOMAINS } from "$projectDir/interproscan/modules/output/representative_domains/main"
include { AGGREGATE_PARSED_SEQS } from "$projectDir/interproscan/modules/output/aggregate_parsed_seqs/main"
include { AGGREGATE_RESULTS } from "$projectDir/interproscan/subworkflows/aggregate_results/main"
include { WRITE_RESULTS } from "$projectDir/interproscan/modules/output/write_results/main"

include { PRE_CHECKS } from "$projectDir/interproscan/subworkflows/pre_checks/main"
include { SEQUENCE_PRECALC } from "$projectDir/interproscan/subworkflows/sequence_precalc/main"
include { SEQUENCE_ANALYSIS } from "$projectDir/interproscan/subworkflows/sequence_analysis/main"
include { XREFS } from "$projectDir/interproscan/subworkflows/xrefs/main"

workflow {
    println "# ${workflow.manifest.name} ${workflow.manifest.version}"
    println "# ${workflow.manifest.description}\n"

    // Help text
    if (params.keySet().any { it.equalsIgnoreCase("help") }) {
        println InterProScan.getHelpMessage(params.appsConfig)
        exit 0
    }

    // Params validation
    InterProScan.checkParams(params)

    // Applications validation
    params.appsToRun = InterProScan.checkApplications(params)

    // Perform preliminary validation checks before running the analysis
    if (params.input != null) {
        input_file = file(params.input)
    } else {
        input_file = null
    }

    /*
    The data dir path is reconfigured and passed from PRE_CHECKS
    to SEQUENCE_ANALYSIS to ensure PRE_CHECKS is completed before
    the other subworkflow starts.
    */
    dataDirPath = Channel.empty()
    PRE_CHECKS(
        params.help,
        input_file,
        params.datadir,
        params.nucleic,
        params.keySet(),
        params.applications,
        params.formats,
        params.version,
        params.ipscn_version,
        params.signalp_mode,
        params.signalp_gpu,
        params.goterms,
        params.pathways
    )
    dataDirPath = file(PRE_CHECKS.out.dataDir.val).toAbsolutePath().toString()
    log.info "Using data files located in ${dataDirPath}"

    applications = (params.applications.toLowerCase().split(',') as Set).join(',')

    Channel.fromPath( input_file , checkIfExists: true)
    .unique()
    .splitFasta( by: params.batchsize, file: true )
    .set { ch_fasta }

    // if nucleic acid seqs provided, predict ORFs
    // either way, then break up input FASTA into batches
    if (params.nucleic) {
        if (params.translate.strand.toLowerCase() !in ['both','plus','minus']) {
            log.info "Strand option '${params.translate.strand.toLowerCase()}' in nextflow.config not recognised. Accepted: 'both', 'plus', 'minus'"
            exit 1
        }
        GET_ORFS(
            ch_fasta,
            params.translate.strand,
            params.translate.methionine,
            params.translate.min_len,
            params.translate.genetic_code
        )
        GET_ORFS.out.splitFasta( by: params.batchsize, file: true )
        .set { orfs_fasta }
        /* Provide the translated ORFs and the original nts seqs
        So that the ORFs can be associated with the source nucleic seq
        in the final output */
        PARSE_SEQUENCE(orfs_fasta, ch_fasta, params.nucleic)
    }
    else {
        PARSE_SEQUENCE(ch_fasta, ch_fasta, params.nucleic)
    }

    disable_precalc = params.disable_precalc
    sequences_to_analyse = null
    parsed_matches = Channel.empty()
    if (!disable_precalc) {
        log.info "Using precalculated match lookup service"
        SEQUENCE_PRECALC(PARSE_SEQUENCE.out, applications, false)  // final: bool to indicate not a unit test
        sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
        parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
    }

    if (parsed_matches.collect() == null) {
            // cases in which the lookup check ran successfully but lookup matches not
            disable_precalc = true
            log.info "ERROR: unable to connect to match lookup service. Max retries reached. Running analysis locally..."
    }

    analysis_result = Channel.empty()
    if (disable_precalc || sequences_to_analyse) {
        log.info "Running sequence analysis"
        if (sequences_to_analyse && !disable_precalc) {
            fasta_to_runner = sequences_to_analyse
        }
        else {
            if (params.nucleic) {
                fasta_to_runner = orfs_fasta
            }
            else {
                fasta_to_runner = ch_fasta
            }
        }
        parsed_analysis = SEQUENCE_ANALYSIS(
            fasta_to_runner,
            applications,
            dataDirPath,
            params.signalp_mode
        )
    }

    AGGREGATE_PARSED_SEQS(PARSE_SEQUENCE.out.collect())

    all_results = parsed_matches.concat(parsed_analysis)
    AGGREGATE_RESULTS(all_results)

    /* XREFS:
    Add signature and entry desc and names
    Add PAINT annotations (if panther is enabled)
    Add go terms (if enabled)
    Add pathways (if enabled)
    */
    XREFS(AGGREGATE_RESULTS.out, applications, dataDirPath)

    REPRESENTATIVE_DOMAINS(XREFS.out.collect())

    Channel.from(params.formats.toLowerCase().split(','))
    .set { ch_format }

    WRITE_RESULTS(
        input_file.getName(),
        AGGREGATE_PARSED_SEQS.out,
        REPRESENTATIVE_DOMAINS.out.collect(),
        ch_format,
        params.outdir,
        params.ipscn_version,
        params.nucleic
    )
}

workflow.onComplete = {
    def input_file = file(params.input)
    def outputFileName = input_file.getName()
    def outputDir = params.outdir.endsWith('/') ? params.outdir[0..-2] : params.outdir

    println "InterProScan workflow completed successfully: $workflow.success."
    println "Any results are located at ${outputDir}/${outputFileName}.ips6.*"
    println "Duration: $workflow.duration"
}
