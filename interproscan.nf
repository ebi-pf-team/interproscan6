nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES AND SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PARSE_SEQUENCE } from "$projectDir/modules/local/parse_sequence/main"
include { GET_ORFS } from "$projectDir/modules/local/get_orfs/main"
include { AGGREGATE_RESULTS } from "$projectDir/modules/local/write_output/aggregate_results/main"

include { PRE_CHECKS } from "$projectDir/subworkflows/pre_checks/main"
include { SEQUENCE_PRECALC } from "$projectDir/subworkflows/sequence_precalc/main"
include { SEQUENCE_ANALYSIS } from "$projectDir/subworkflows/sequence_analysis/main"


workflow {
    // Perform preliminary validation checks before running the analysis
    PRE_CHECKS(params)
    System.exit(0)

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
