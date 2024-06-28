include {
    CDD_RUNNER;
    CDD_POSTPROCESS;
    CDD_PARSER
} from "$projectDir/modules/cdd/main"
include {
    FUNFAM_HMMER_RUNNER;
    HMMER_RUNNER as GENERIC_HMMER_RUNNER;
    HMMER_RUNNER as GENE3D_HMMER_RUNNER;
    HMMER_RUNNER as HAMAP_HMMER_RUNNER;
    HMMER_RUNNER as SFLD_HMMER_RUNNER;
    HMMER_RUNNER as PANTHER_HMMER_RUNNER;
} from "$projectDir/modules/hmmer/runner/main"
include {
    HMMER_PARSER as GENERIC_HMMER_PARSER;
    HMMER_PARSER as FUNFAM_HMMER_PARSER;
    HMMER_PARSER as GENE3D_HMMER_PARSER;
    HMMER_PARSER as HAMAP_HMMER_PARSER;
    HMMER_PARSER as SFLD_HMMER_PARSER;
    HMMER_PARSER as PANTHER_HMMER_PARSER;
} from "$projectDir/modules/hmmer/parser/main"
include {
    CATH_RESEOLVE_HITS as FUNFAM_CATH_RESEOLVE_HITS;  // third party tool to minimise suprious hits
    ADD_CATH_SUPERFAMILIES as FUNFAM_ADD_CATH_SUPERFAMILIES;  // used for gene3D and Funfam
    CATH_RESEOLVE_HITS as GENE3D_CATH_RESEOLVE_HITS;
    ADD_CATH_SUPERFAMILIES as GENE3D_ADD_CATH_SUPERFAMILIES;
    HAMAP_POST_PROCESSER;
    PANTHER_POST_PROCESSER;
    SFLD_POST_PROCESSER
} from "$projectDir/modules/hmmer/post_processing/main"
include {
    FUNFAM_FILTER_MATCHES;
    GENE3D_FILTER_MATCHES;
    HAMAP_FILTER_MATCHES;
    PANTHER_FILTER_MATCHES;
    SFLD_FILTER_MATCHES;
} from "$projectDir/modules/hmmer/filter/main"
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
    boolean gene3d_funfam_processed = false
    // To prevent duplication if Gene3D and Funfam are called

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
                false,  // don't build an alignment file
                false,   // don't build a hmmer.tbl file path 
                []  // no post-processing params
            ]

        /*
        Place FunFam inside the Gene3D post-processing
        because it must run after the Gene3D path
        */
        gene3d_funfam: (member == 'gene3d' || member == 'funfam') && !gene3d_funfam_processed
            gene3d_funfam_processed = true
            return [
                params.members."gene3d".hmm,
                params.members."gene3d".switches,
                params.members."gene3d".release,
                true,   // build an alignment file
                false,   // don't build a hmmer.tbl file path
                [
                    params.members."gene3d".postprocess.cath_resolve_hits_switches,
                    params.members."gene3d".postprocess.model2sf_map,
                    params.members."gene3d".postprocess.discontinuous_regs,
                    params.members."gene3d".postprocess.assign_cath_superfamilies,
                    params.members."funfam".hmm,
                    params.members."funfam".switches,
                    params.members."funfam".release,
                ]
            ]
        
        hamap: runner == 'hamap'
            return [
                params.members."${member}".hmm,
                params.members."${member}".switches,
                params.members."${member}".release,
                false,  // don't build an alignment file
                true,    // build a hmmer.tbl file path
                [
                    params.members."${member}".postprocess.models_dir,
                    params.members."${member}".postprocess.pfsearchv3_switches,
                ]
            ]

        panther: runner == 'panther'
            return [
                params.members."${member}".hmm,
                params.members."${member}".switches,
                params.members."${member}".release,
                false,  // don't build an alignment file
                false,   // don't build a hmmer.tbl file path
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
                true,  // build an alignment file
                false,  // don't build a hmmer.tbl file path
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
                    params.members."${member}".postprocess.switches,
                    params.members."${member}".postprocess.data,
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
    }.set { member_params }

    /*
    Member databases that use HMMER
    */

    // AntiFam and NCBIfam
    runner_hmmer_params = fasta.combine(member_params.hmmer)
    GENERIC_HMMER_RUNNER(runner_hmmer_params)
    GENERIC_HMMER_PARSER(
        GENERIC_HMMER_RUNNER.out[0],  // hmmer.out path
        GENERIC_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        GENERIC_HMMER_RUNNER.out[2],  // post-processing-params
        tsv_pro,
        "false"
    )

    // Cath-Gene3D (+ cath-resolve-hits + assing-cath-superfamilies)
    // These also run for FunFam as Gene3D must be run before FunFam
    runner_gene3d_params = fasta.combine(member_params.gene3d_funfam)
    GENE3D_HMMER_RUNNER(runner_gene3d_params)
    GENE3D_HMMER_PARSER(
        GENE3D_HMMER_RUNNER.out[0],  // hmmer.out path
        GENE3D_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        GENE3D_HMMER_RUNNER.out[2],  // post-processing-params
        tsv_pro,
        "false"
    )
    GENE3D_CATH_RESEOLVE_HITS(
        GENE3D_HMMER_RUNNER.out[0],  // hmmer.out path
        GENE3D_HMMER_RUNNER.out[2],  // post-processing-params
    )
    GENE3D_ADD_CATH_SUPERFAMILIES(GENE3D_CATH_RESEOLVE_HITS.out)
    GENE3D_FILTER_MATCHES(GENE3D_ADD_CATH_SUPERFAMILIES.out, GENE3D_HMMER_PARSER.out)

    // FunFam (+ gene3D + cath-resolve-hits + assing-cath-superfamilies)
    // split into a channel so Nextflow can automatically manage the parallel execution of HmmSearch
    GENE3D_FILTER_MATCHES.out[1]
        .splitText() { it.replace('\n', '') }
        .set { funfam_cath_superfamilies }
    runner_funfam_params = fasta.combine(member_params.gene3d_funfam)
    runner_funfam_params_with_cath = runner_funfam_params.combine(funfam_cath_superfamilies)

    FUNFAM_HMMER_RUNNER(runner_funfam_params_with_cath, applications)
    FUNFAM_HMMER_PARSER(FUNFAM_HMMER_RUNNER.out, tsv_pro, "false")
    FUNFAM_CATH_RESEOLVE_HITS(
        FUNFAM_HMMER_RUNNER.out[0],  // hmmer.out path
        FUNFAM_HMMER_RUNNER.out[2]   // post-processing-params
    )
    FUNFAM_FILTER_MATCHES(FUNFAM_HMMER_PARSER.out, FUNFAM_CATH_RESEOLVE_HITS.out)

    // HAMAP (+ pfsearch_wrapper.py)
    runner_hamap_params = fasta.combine(member_params.hamap)
    HAMAP_HMMER_RUNNER(runner_hamap_params)
    HAMAP_HMMER_PARSER(
        HAMAP_HMMER_RUNNER.out[0],  // hmmer.out path
        HAMAP_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        HAMAP_HMMER_RUNNER.out[2],  // post-processing-params
        tsv_pro,
        "false"
    )
    HAMAP_POST_PROCESSER(
        HAMAP_HMMER_RUNNER.out[5], // path to fasta file parsed by HMMER
        HAMAP_HMMER_RUNNER.out[4], // hmmer .tbl file path
        HAMAP_HMMER_RUNNER.out[2]  // post-processing-params
    )
    HAMAP_FILTER_MATCHES(
        HAMAP_HMMER_PARSER.out,    // internal IPS6 JSON
        HAMAP_POST_PROCESSER.out   // output from pfsearch_wrapper.py
    )

    // Panther (+ treegrafter + epa-ng)
    runner_panther_params = fasta.combine(member_params.panther)
    PANTHER_HMMER_RUNNER(runner_panther_params)
    PANTHER_HMMER_PARSER(
        PANTHER_HMMER_RUNNER.out[0],  // hmmer.out path
        PANTHER_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        PANTHER_HMMER_RUNNER.out[2],  // post-processing-params
        tsv_pro,
        "false"
    )
    PANTHER_POST_PROCESSER(
        PANTHER_HMMER_PARSER.out,     // hmmer.out path
        PANTHER_HMMER_RUNNER.out[2],  // post-processing-params
        fasta
    )
    PANTHER_FILTER_MATCHES(
        PANTHER_HMMER_PARSER.out,   // internal ips6 json
        PANTHER_POST_PROCESSER.out  // treegrafter output + post-processing params
    )

    // SFLD (+ post-processing binary to add sites and filter hits)
    runner_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_sfld_params)
    SFLD_HMMER_PARSER(
        SFLD_HMMER_RUNNER.out[0],  // hmmer.out path
        SFLD_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        SFLD_HMMER_RUNNER.out[2],  // post-processing-params
        tsv_pro,
        "true"
    )
    SFLD_POST_PROCESSER(
        SFLD_HMMER_RUNNER.out[0],  // hmmer.out path
        SFLD_HMMER_RUNNER.out[1],  // hmmer.dtbl path
        SFLD_HMMER_RUNNER.out[2],  // post-processing-params
        SFLD_HMMER_RUNNER.out[3],  // alignment file
        tsv_pro,
    )
    SFLD_FILTER_MATCHES(SFLD_HMMER_PARSER.out, SFLD_POST_PROCESSER.out)

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

    if (applications.contains("gene3d")) {
        GENERIC_HMMER_PARSER.out[0].concat(
            FUNFAM_FILTER_MATCHES.out[0],
            GENE3D_FILTER_MATCHES.out[0],
            HAMAP_FILTER_MATCHES.out,
            PANTHER_FILTER_MATCHES.out,
            SFLD_FILTER_MATCHES.out,
            CDD_PARSER.out,
            SIGNALP_PARSER.out
        )
        .set { parsed_results }
    }
    else {
        GENERIC_HMMER_PARSER.out[0].concat(
            FUNFAM_FILTER_MATCHES.out[0],
            HAMAP_FILTER_MATCHES.out,
            PANTHER_FILTER_MATCHES.out,
            SFLD_FILTER_MATCHES.out,
            CDD_PARSER.out,
            SIGNALP_PARSER.out
        )
        .set { parsed_results }
    }

    emit:
    parsed_results
}
