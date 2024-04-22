include { HMMER_RUNNER as GENERIC_HMMER_RUNNER; HMMER_RUNNER as SFLD_HMMER_RUNNER } from "$projectDir/modules/local/hmmer/runner/main"
include { HMMER_PARSER } from "$projectDir/modules/local/hmmer/parser/main"
include { GENE3D_POST_PROCESSER; FUNFAM_POST_PROCESSER; SFLD_POST_PROCESSER } from "$projectDir/modules/local/hmmer/post_processing/main"
include { SFLD_PARSER } from "$projectDir/modules/local/hmmer/parser/slfd"
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
        // funfam
        // gene3d
        // hamap
        // panther
        // pfam
        // pirsf
        // pirsr
        // smart ?
        // superfamily
        if (member == 'sfld') {
            runner = 'sfld'
        }

        /*
        The post processing of some applications (e.g. SFLD) hits requires additional files
        and parameters relative to the generic hmmer runner and parser
        */
        hmmer: runner == 'hmmer'
            return [ 
                params.members."${member}".hmm, params.members."${member}".switches, 
                false, []
            ]
            
        sfld: runner == 'sfld'
            return [
                params.members."${member}".hmm, params.members."${member}".switches,
                true, [
                    params.members."${member}".postprocess.bin,
                    params.members."${member}".postprocess.sites_annotation,
                    params.members."${member}".postprocess.hierarchy
                ]
            ]
        if (params.members."${member}".runner == "signalp") {
            runner = 'signalp'
        }

        signalp: runner == 'signalp'
            return [
                params.members.signalp.data.mode,
                params.members.signalp.data.model_dir,
                params.members.signalp.data.organism,
                params.members.signalp.switches,
                params.members.signalp.data.pvalue
            ]

        other: true
            log.info "Application ${member} (still) not supported"
            
        log.info "Running $runner for $member"
    }.set { member_params }

    runner_hmmer_params = fasta.combine(member_params.hmmer)
    GENERIC_HMMER_RUNNER(runner_hmmer_params)
    HMMER_PARSER(GENERIC_HMMER_RUNNER.out, params.tsv_pro)

    runner_hmmer_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_hmmer_sfld_params)
    SFLD_POST_PROCESSER(SFLD_HMMER_RUNNER.out, params.tsv_pro)
    SFLD_PARSER(SFLD_POST_PROCESSER.out, params.tsv_pro, true)  // set sites to true for SFLD

    runner_signalp_params = fasta.combine(member_params.signalp)
    SIGNALP_RUNNER(runner_signalp_params)
    SIGNALP_PARSER(SIGNALP_RUNNER.out, params.tsv_pro)

    HMMER_PARSER.out.concat(SIGNALP_PARSER.out)
    .set { parsed_results }

    emit:
    parsed_results
}
