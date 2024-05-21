include {
    CDD_RUNNER;
    CDD_POSTPROCESS;
    CDD_PARSER
} from "$projectDir/modules/cdd/main"
include {
    HMMER_RUNNER as GENERIC_HMMER_RUNNER;
    HMMER_RUNNER as SFLD_HMMER_RUNNER;
    HMMER_RUNNER as PANTHER_HMMER_RUNNER;
} from "$projectDir/modules/hmmer/runner/main"
include {
    HMMER_PARSER as GENERIC_HMMER_PARSER;
    HMMER_PARSER as SFLD_HMMER_PARSER;
    HMMER_PARSER as PANTHER_HMMER_PARSER;
} from "$projectDir/modules/hmmer/parser/main"
include {
    PANTHER_POST_PROCESSER;
    SFLD_POST_PROCESSER
} from "$projectDir/modules/hmmer/post_processing/main"
include {
    SIGNALP_RUNNER;
    SIGNALP_PARSER
 } from "$projectDir/modules/signalp/main"


workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications
    tsv_pro

    main:
    // Divide members up into their respective analysis pipelines/methods
    Channel.from(applications.split(','))
    .branch { member ->
        release = params.members."${member}".release
        log.info "Running $member version $release"
        runner = ''

        if (member == 'antifam' || member == "ncbifam") {
            runner = 'hmmer'
        } else {
            runner = member
        }

        /*
        Member databases that use HMMER:
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
                    params.members."${member}".postprocess.paint_annotations,
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

        /*
        Member databases that do NOT use HMMER
        */

        cdd: runner == "cdd"
            return [
                params.members."${member}".library,
                params.members."${member}".release,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.bin,
                    params.members."${member}".postprocess.switches,
                    params.members."${member}".postprocess.data,
                    params.members."${member}".postprocess.signature_list
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
    }.set { member_params }

    /*
    Member databases that use HMMER
    */

    // AntiFam and NCBIfam
    runner_hmmer_params = fasta.combine(member_params.hmmer)
    GENERIC_HMMER_RUNNER(runner_hmmer_params)
    GENERIC_HMMER_PARSER(GENERIC_HMMER_RUNNER.out, tsv_pro, "antifam")

    // Panther (+ treegrafter + epa-ng)
    runner_hmmer_panther_params = fasta.combine(member_params.panther)
    PANTHER_HMMER_RUNNER(runner_hmmer_panther_params)
    PANTHER_HMMER_PARSER(PANTHER_HMMER_RUNNER.out, tsv_pro, "panther")
    PANTHER_POST_PROCESSER(PANTHER_HMMER_PARSER.out, fasta)

    // SFLD (+ post-processing binary to add sites and filter hits)
    runner_hmmer_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_hmmer_sfld_params)
    SFLD_HMMER_PARSER(SFLD_HMMER_RUNNER.out, tsv_pro, "sfld")
    SFLD_POST_PROCESSER(SFLD_HMMER_PARSER.out, tsv_pro)
    
    /*
    Member databases that do NOT use HMMER
    */

    // CDD
    runner_cdd_params = fasta.combine(member_params.cdd)
    CDD_RUNNER(runner_cdd_params)
    CDD_POSTPROCESS(CDD_RUNNER.out)
    CDD_PARSER(CDD_POSTPROCESS.out)

    // SignalP
    runner_signalp_params = fasta.combine(member_params.signalp)
    SIGNALP_RUNNER(runner_signalp_params)
    SIGNALP_PARSER(SIGNALP_RUNNER.out, tsv_pro)

    /*
    Gather the results
    */

    GENERIC_HMMER_PARSER.out.concat(
        PANTHER_HMMER_PARSER.out,
        SFLD_HMMER_PARSER.out,
        CDD_PARSER.out,
        SIGNALP_PARSER.out
    )
    .set { parsed_results }  // gathers the paths of the output file from each process

    emit:
    parsed_results
}
