include { 
    HMMER_RUNNER as GENERIC_HMMER_RUNNER;
    HMMER_RUNNER as SFLD_HMMER_RUNNER;
    HMMER_RUNNER as PANTHER_HMMER_RUNNER; 
} from "$projectDir/modules/local/hmmer/runner/main"
include { 
    HMMER_PARSER as GENERIC_HMMER_PARSER;
    HMMER_PARSER as SFLD_HMMER_PARSER;
    HMMER_PARSER as PANTHER_HMMER_PARSER;
} from "$projectDir/modules/local/hmmer/parser/main"
include { 
    PANTHER_POST_PROCESSER;
    SFLD_POST_PROCESSER
} from "$projectDir/modules/local/hmmer/post_processing/main"
include { SIGNALP_RUNNER } from "$projectDir/modules/local/signalp/runner/main"
include { SIGNALP_PARSER } from "$projectDir/modules/local/signalp/parser/main"


workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications

    main:
    // Divide members up into their respective analysis pipelines/methods
    Channel.from(applications.split(','))
    .branch { member ->
        runner = ''
        if (params.members."${member}".runner == "hmmer") {
            runner = 'hmmer'
        }
        if (member == 'panther') {
            runner = 'panther'
        } else if (member == 'sfld') {
            runner = 'sfld'
        } else if (member == 'signalp') {
            runner = 'signalp'
        }

        /*
        The post processing of some applications (e.g. SFLD) hits requires additional files
        and parameters relative to the generic hmmer runner and parser
        */
        hmmer: runner == 'hmmer'
            return [
                params.members."${member}".hmm,
                params.members."${member}".switches,
                params.members."${member}".release,
                false, []
            ]

        panther: runner == 'panther'
            return [
                params.members."${member}".hmm,
                params.members."${member}".switches,
                params.members."${member}".release,
                false,
                [
                    params.members."${member}".postprocess.data_dir,
                    params.members."${member}".postprocess.evalue,
                ]
            ]

        sfld: runner == 'sfld'
            return [
                params.members."${member}".hmm,
                params.members."${member}".switches,
                params.members."${member}".release,
                true,
                [
                    params.members."${member}".postprocess.bin,
                    params.members."${member}".postprocess.sites_annotation,
                    params.members."${member}".postprocess.hierarchy
                ]
            ]

        signalp: runner == 'signalp'
            return [
                params.members.signalp.data.mode,
                params.members.signalp.data.model_dir,
                params.members.signalp.data.organism,
                params.members.signalp.switches,
                params.members.signalp.data.pvalue,
                params.members.signalp.release
            ]

        other: true
            log.info "Application ${member} (still) not supported"

        log.info "Running $runner for $member"
    }.set { member_params }

    // AntiFam and NCBIfam
    runner_hmmer_params = fasta.combine(member_params.hmmer)
    GENERIC_HMMER_RUNNER(runner_hmmer_params)
    GENERIC_HMMER_PARSER(GENERIC_HMMER_RUNNER.out, params.tsv_pro, false)  // set sites to false

    // Panther (+ treegrafter + epa-ng)
    runner_hmmer_panther_params = fasta.combine(member_params.panther)
    PANTHER_HMMER_RUNNER(runner_hmmer_panther_params)
    PANTHER_POST_PROCESSER(PANTHER_HMMER_RUNNER.out, fasta)
    PANTHER_HMMER_PARSER(PANTHER_POST_PROCESSER.out, params.tsv_pro, false)

    // SFLD (+ post-processing binary to add sites and filter hits)
    runner_hmmer_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_hmmer_sfld_params)
    SFLD_POST_PROCESSER(SFLD_HMMER_RUNNER.out, params.tsv_pro)
    SFLD_HMMER_PARSER(SFLD_POST_PROCESSER.out, params.tsv_pro, true)  // set sites to true for SFLD

    // SignalP
    runner_signalp_params = fasta.combine(member_params.signalp)
    SIGNALP_RUNNER(runner_signalp_params)
    SIGNALP_PARSER(SIGNALP_RUNNER.out, params.tsv_pro)

    GENERIC_HMMER_PARSER.out.concat(
        PANTHER_HMMER_PARSER.out,
        SFLD_HMMER_PARSER.out,
        SIGNALP_PARSER.out
    )
    .set { parsed_results }  // gathers the paths of the output file from each process

    emit:
    parsed_results
}
