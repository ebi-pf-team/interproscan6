nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/parse_sequence/main"
include { GET_ORFS } from "$projectDir/modules/get_orfs/main"
include { AGGREGATE_RESULTS } from "$projectDir/modules/output/aggregate_results/main"
include { WRITE_RESULTS } from "$projectDir/modules/output/write_results/main"

include { PRE_CHECKS } from "$projectDir/subworkflows/pre_checks/main"
include { SEQUENCE_PRECALC } from "$projectDir/subworkflows/sequence_precalc/main"
include { SEQUENCE_ANALYSIS } from "$projectDir/subworkflows/sequence_analysis/main"
include { XREFS } from "$projectDir/subworkflows/xrefs/main"

workflow {
    // Perform preliminary validation checks before running the analysis
    PRE_CHECKS(
        params.help,
        file(params.input),
        params.nucleic,
        params.keySet(),
        params.applications,
        params.formats
    )

    formats = params.formats.toLowerCase()

    tsv_pro = false
    if (formats.contains("tsv-pro")) {
        tsv_pro = true
    }

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
        SEQUENCE_PRECALC(PARSE_SEQUENCE.out, applications, false)  // final: bool to indicate not a unit test
        sequences_to_analyse = SEQUENCE_PRECALC.out.sequences_to_analyse
        parsed_matches = SEQUENCE_PRECALC.out.parsed_matches
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
        parsed_analysis = SEQUENCE_ANALYSIS(fasta_to_runner, applications, tsv_pro)
    }

    all_results = parsed_matches.concat(parsed_analysis)

    AGGREGATE_RESULTS(all_results.collect())

    XREFS(AGGREGATE_RESULTS.out)

    Channel.from(formats.split(','))
    .set { ch_format }

    WRITE_RESULTS(PARSE_SEQUENCE.out.collect(), XREFS.out.collect(), ch_format, params.output, params.ipsc_version)
}

workflow.onComplete = {
    println "InterProScan workflow completed $workflow.success. Results in ${params.output}.*"
    println "Duration: $workflow.duration"
}

log.info """
If you use InterProScan in your work please cite:

InterProScan:
> Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H,
Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A,
Scheremetjew M, Yong SY, Lopez R, Hunter S.
InterProScan 5: genome-scale protein function classification.
Bioinformatics. 2014 May 1;30(9):1236-40. doi: 10.1093/bioinformatics/btu031.
Epub 2014 Jan 21. PMID: 24451626; PMCID: PMC3998142.

InterPro:
> Paysan-Lafosse T, Blum M, Chuguransky S, Grego T, Pinto BL, Salazar GA, Bileschi ML,
Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunić I, Marchler-Bauer A, Mi H,
Natale DA, Orengo CA, Pandurangan AP, Rivoire C, Sigrist CJA, Sillitoe I, Thanki N,
Thomas PD, Tosatto SCE, Wu CH, Bateman A.
InterPro in 2022. Nucleic Acids Res. 2023 Jan 6;51(D1):D418-D427.
doi: 10.1093/nar/gkac993. PMID: 36350672; PMCID: PMC9825450.
"""
