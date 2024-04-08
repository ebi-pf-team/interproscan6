include { HMMER_RUNNER } from "$projectDir/modules/local/hmmer/runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer/parser/main"
include { SIGNALP_RUNNER } from "$projectDir/modules/local/signalp/runner/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    log.info "DEBUG params: ${params.members.signalp.data.mode},${params.members.signalp.data.model_dir},${params.members.signalp.data.organism},${params.members.signalp.switches}"
    Channel.from(applications.split(','))
    .branch { member ->
        def runner = ''
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }
        if (params.members."${member}".runner == "signalp") {
            runner = 'signalp'
        }

        log.info "Running ${runner}"
        hmmer: runner == 'hmmer'
            return [
                params.members."${member}".data,
                params.members."${member}".switches
            ]
        signalp: runner == 'signalp'
            return [
                params.members.signalp.data.mode,
                params.members.signalp.data.model_dir,
                params.members.signalp.data.organism,
                params.members.signalp.switches
            ]
        other: true
            log.info "Application ${member} (still) not supported"
    }.set { member_params }

    if (runner == 'hmmer') {
        runner_hmmer_params = fasta.combine(member_params.hmmer)
        HMMER_RUNNER(runner_hmmer_params)
        HMMER_PARSER(HMMER_RUNNER.out, params.tsv_pro)
        result = HMMER_PARSER.out
    }
    if (runner == 'signalp') {
        runner_signalp_params = fasta.combine(member_params.signalp)
        log.info "signalp params: ${runner_signalp_params}"
        SIGNALP_RUNNER(runner_signalp_params)
//         SIGNALP_PARSER(...)
//         result = SIGNALP_PARSER.out
        result = SIGNALP_RUNNER.out
    }

    emit:
    result
}
