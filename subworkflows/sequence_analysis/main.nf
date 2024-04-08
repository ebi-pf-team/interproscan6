include { HMMER_RUNNER } from "$projectDir/modules/local/hmmer/runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer/parser/main"
include { SIGNALP_RUNNER } from "$projectDir/modules/local/signalp/runner/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    Channel.from(applications.split(','))
    .branch { member ->
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }
        if (params.members."${member}".runner == "signalp") {
            runner = 'signalp'
        }

        log.info "runner: $runner"
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

<<<<<<< HEAD
    log.info "Running HMMER"
    runner_hmmer_params = fasta.combine(member_params.hmmer)
    HMMER_RUNNER(runner_hmmer_params)
    HMMER_PARSER(HMMER_RUNNER.out, params.tsv_pro)

    log.info "Running SignalP"
    runner_signalp_params = fasta.combine(member_params.signalp)
    SIGNALP_RUNNER(runner_signalp_params)

    HMMER_PARSER.out.concat(SIGNALP_RUNNER.out)
    .set { parsed_results }
=======
    if (runner == 'hmmer') {
        runner_hmmer_params = fasta.combine(member_params.hmmer)
        HMMER_RUNNER(runner_hmmer_params)
        HMMER_PARSER(HMMER_RUNNER.out, params.tsv_pro)
        return HMMER_PARSER.out
    }.set { result }
    else
    if (runner == 'signalp') {
        runner_signalp_params = fasta.combine(member_params.signalp)
        SIGNALP_RUNNER(runner_signalp_params)
//         SIGNALP_PARSER(...)
//         result = SIGNALP_PARSER.out
        return SIGNALP_RUNNER.out
    }.set { result }
    else {
        log.info "Runner ${runner} (still) not supported"
    }
>>>>>>> 6772592 (removing logs)

    emit:
    parsed_results
}
