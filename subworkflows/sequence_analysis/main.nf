include { HMMER_RUNNER } from "$projectDir/modules/local/hmmer_runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer_parser/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    Channel.from(applications.split(','))
    .branch { member ->
        def runner = ''
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }

        hmmer: runner == 'hmmer'
            return [ params.members."${member}".data, params.members."${member}".switches ]
        other: true
            log.info "Application ${member} (still) not supported"
    }.set { member_params }

    runner_hmmer_params = fasta.combine(member_params.hmmer)
    HMMER_RUNNER(runner_hmmer_params)
    HMMER_PARSER(HMMER_RUNNER.out, params.tsv_pro)

    emit:
    HMMER_PARSER.out
}
