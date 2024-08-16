include {
    CDD_RUNNER;
    CDD_PARSER;
    CDD_POSTPROCESS;
} from "$projectDir/interproscan/modules/cdd/main"
include {
    COILS_RUNNER;
    COILS_PARSER;
} from "$projectDir/interproscan/modules/coils/main"
include {
    HMMER_RUNNER as ANTIFAM_HMMER_RUNNER;
    HMMER_RUNNER as NCBIFAM_HMMER_RUNNER;
    HMMER_RUNNER as GENE3D_HMMER_RUNNER;
    HMMER_RUNNER as PANTHER_HMMER_RUNNER;
    HMMER_RUNNER as PFAM_HMMER_RUNNER;
    HMMER_RUNNER as PIRSR_HMMER_RUNNER;
    HMMER_RUNNER_WITH_ALIGNMENTS as SFLD_HMMER_RUNNER;
    HMMER_SCAN_RUNNER as SUPERFAMILY_HMMER_RUNNER;
    FUNFAM_HMMER_RUNNER;
    HAMAP_HMMER_RUNNER;
    PIRSF_HMMER_RUNNER;
    SMART_HMMER2_RUNNER;
    HMMER_SCAN_RUNNER;
} from "$projectDir/interproscan/modules/hmmer/runner/main"
include {
    HMMER_PARSER as ANTIFAM_HMMER_PARSER;
    HMMER_PARSER as GENE3D_HMMER_PARSER;
    HMMER_PARSER as HAMAP_HMMER_PARSER;
    HMMER_PARSER as NCBIFAM_HMMER_PARSER;
    HMMER_PARSER as PANTHER_HMMER_PARSER;
    HMMER_PARSER as PFAM_HMMER_PARSER;
    HMMER_PARSER as PIRSR_HMMER_PARSER;
    HMMER_PARSER as SFLD_HMMER_PARSER;
    FUNFAM_HMMER_PARSER;
    HMMER_SCAN_PARSER as PIRSF_HMMER_PARSER;
    HMMER2_PARSER;
} from "$projectDir/interproscan/modules/hmmer/parser/main"
include {
    CATH_RESOLVE_HITS as GENE3D_CATH_RESOLVE_HITS;
    ADD_CATH_SUPERFAMILIES as GENE3D_ADD_CATH_SUPERFAMILIES;
    FUNFAM_CATH_RESOLVE_HITS;
    HAMAP_POST_PROCESSER;
    PANTHER_POST_PROCESSER;
    SFLD_POST_PROCESSER;
    SUPERFAMILY_POST_PROCESSER;
} from "$projectDir/interproscan/modules/hmmer/post_processing/main"
include {
    FUNFAM_FILTER_MATCHES;
    GENE3D_FILTER_MATCHES;
    HAMAP_FILTER_MATCHES;
    PANTHER_FILTER_MATCHES;
    PFAM_FILTER_MATCHES;
    PIRSF_FILTER_MATCHES;
    PIRSR_FILTER_MATCHES;
    SFLD_FILTER_MATCHES;
    SMART_FILTER_MATCHES;
    SUPERFAMILY_FILTER_MATCHES;
} from "$projectDir/interproscan/modules/hmmer/filter/main"
include {
    MOBIDB_RUNNER;
    MOBIDB_PARSER;
} from "$projectDir/interproscan/modules/mobidb/main"
include {
    PRINTS_RUNNER;
    PRINTS_PARSER;
} from "$projectDir/interproscan/modules/prints/main"
include {
    PHOBIUS_RUNNER;
    PHOBIUS_PARSER;
} from "$projectDir/interproscan/modules/phobius/main"
include {
    PFSEARCH_RUNNER as PROSITE_PROFILES_RUNNER
} from "$projectDir/interproscan/modules/prosite/pfsearch/runner/main"
include {
    PFSEARCH_PARSER as PROSITE_PROFILES_PARSER
} from "$projectDir/interproscan/modules/prosite/pfsearch/parser/main"
include {
    PFSCAN_RUNNER as PROSITE_PATTERNS_RUNNER
} from "$projectDir/interproscan/modules/prosite/pfscan/runner/main"
include {
    PFSCAN_PARSER as PROSITE_PATTERNS_PARSER
} from "$projectDir/interproscan/modules/prosite/pfscan/parser/main"
include {
    SIGNALP_RUNNER;
    SIGNALP_PARSER;
    SIGNALP_RUNNER as SIGNALP_EUK_RUNNER;
    SIGNALP_PARSER as SIGNALP_EUK_PARSER;
} from "$projectDir/interproscan/modules/signalp/main"

workflow SEQUENCE_ANALYSIS {
    take:
    fasta
    applications
    signalp_mode

    main:
    boolean gene3d_funfam_processed = false
    // To prevent duplication if Gene3D and Funfam are called

    // Divide members up into their respective analysis pipelines/methods
    Channel.from(applications.split(','))
    .branch { member ->
        runner = params.members."${member}".runner
//         log.info "Running $member version $release"

        /*
        Member databases that use HMMER:
        The post processing of some applications (e.g. SFLD) hits requires additional files
        and parameters relative to the generic hmmer runner and parser
        */
        antifam: member == 'antifam'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [] // no post-processing
            ]

        /*
        Place FunFam inside the Gene3D post-processing
        because it must run after the Gene3D path
        */
        gene3d_funfam: (member == 'gene3d' || member == 'funfam') && !gene3d_funfam_processed
            gene3d_funfam_processed = true
            return [
                "gene3d",
                params.members."gene3d".hmm,
                params.members."gene3d".switches,
                [
                    params.members."gene3d".postprocess.cath_resolve_hits_switches,
                    params.members."gene3d".postprocess.model2sf_map,
                    params.members."gene3d".postprocess.discontinuous_regs,
                    params.members."gene3d".postprocess.assign_cath_superfamilies,
                    params.members."funfam".hmm,
                    params.members."funfam".switches,
                ]
            ]

        hamap: member == 'hamap'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.models_dir,
                    params.members."${member}".postprocess.pfsearchv3_switches,
                ]
            ]

        ncbifam: member == 'ncbifam'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [] // no post-processing
            ]

        panther: member == 'panther'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.data_dir,
                    params.members."${member}".postprocess.evalue,
                    params.members."${member}".postprocess.paint_annotations,
                ]
            ]

        pfam: member == 'pfam'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.min_length,
                    params.members."${member}".postprocess.seed,
                    params.members."${member}".postprocess.clan,
                    params.members."${member}".postprocess.data,
                ]
            ]

        pirsf: member == 'pirsf'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.data
                ]
            ]

       pirsr: member == 'pirsr'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.rules
                ]
            ]

       sfld: member == 'sfld'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.sites_annotation,
                    params.members."${member}".postprocess.hierarchy,
                ]
            ]

        superfamily: member == 'superfamily'
        return [
            "${member}",
            params.members."${member}".hmm,
            params.members."${member}".switches,
            [
                params.members."${member}".postprocess.bin,
                params.members."${member}".postprocess.self_hits,
                params.members."${member}".postprocess.cla,
                params.members."${member}".postprocess.model,
                params.members."${member}".postprocess.pdbj95d,
                params.members."${member}".postprocess.ass3_switches,
            ]
        ]

        // uses HMMER2, has a slightly different set up
        smart: member == 'smart'
            return [
                "${member}",
                params.members."${member}".hmm,
                params.members."${member}".switches,
            ]

        /*
        Member databases that do NOT use HMMER
        */
        cdd: member == "cdd"
            return [
                params.members."${member}".library,
                params.members."${member}".switches,
                [
                    params.members."${member}".postprocess.switches,
                    params.members."${member}".postprocess.data,
                ]
            ]
        coils: member == "coils"
            return [
            params.members."${member}".switches
            ]

        mobidb: member == "mobidb"
            return [
                params.members."${member}".switches
            ]

        phobius: member == "phobius"
            return []

        prints: member == 'prints'
            return [
                params.members.prints.data.hierarchy,
                params.members.prints.data.pval,
                params.members.prints.switches
            ]

        prosite_patterns: member == "prosite_patterns"
            return [
                params.members."${member}".data,
                params.members."${member}".evaluator,
                params.members."${member}".switches
            ]

        prosite_profiles: member == "prosite_profiles"
            return [
                params.members."${member}".data,
                params.members."${member}".switches,
                params.members."${member}".skip_flagged_profiles
            ]

        signalp: member == 'signalp'

            return [
                params.signalp_mode,
                params.members.signalp.data.organism,
                params.members.signalp.data.pvalue,
            ]

        signalp_euk: member == 'signalp_euk'

            return [
                params.signalp_mode,
                params.members.signalp_euk.data.model_dir,
                params.members.signalp_euk.data.organism,
                params.members.signalp_euk.switches,
                params.members.signalp_euk.data.pvalue,
                params.members.signalp_euk.release
            ]
    }.set { member_params }

    /*
    Member databases that use HMMER
    */

    /*
    Using generic HMMER
    */

    // AntiFam
    runner_antifam_params = fasta.combine(member_params.antifam)
    ANTIFAM_HMMER_RUNNER(runner_antifam_params)
    ANTIFAM_HMMER_PARSER(ANTIFAM_HMMER_RUNNER.out[0])  // hmmer.out path

    // Cath-Gene3D (+ cath-resolve-hits + assing-cath-superfamilies)
    // These also run for FunFam as Gene3D must be run before FunFam
    runner_gene3d_params = fasta.combine(member_params.gene3d_funfam)
    GENE3D_HMMER_RUNNER(runner_gene3d_params)
    GENE3D_HMMER_PARSER(GENE3D_HMMER_RUNNER.out[0])  // hmmer.out path
    GENE3D_CATH_RESOLVE_HITS(GENE3D_HMMER_RUNNER.out)
    GENE3D_ADD_CATH_SUPERFAMILIES(
        GENE3D_CATH_RESOLVE_HITS.out, // cath-resolve-hits out file
        GENE3D_HMMER_RUNNER.out[1]    // post-processing-params
    )
    GENE3D_FILTER_MATCHES(
        GENE3D_ADD_CATH_SUPERFAMILIES.out,  // add-superfams out file
        GENE3D_HMMER_PARSER.out,            // ips6 json
        GENE3D_HMMER_RUNNER.out[1]          // post-processing-params
    )
    // Gene3D filter matches out puts (0) ips6 json; (1) cath superfamilies txt file

    // FunFam (+ gene3D + cath-resolve-hits + assing-cath-superfamilies)
    // split into a channel so Nextflow can automatically manage the parallel execution of HmmSearch
    runner_funfam_params = fasta.combine(member_params.gene3d_funfam)
    FUNFAM_HMMER_RUNNER(
        runner_funfam_params,              // hmmer runner input tuple
        GENE3D_FILTER_MATCHES.out[1],      // cath_superfamilies txt file
        applications                       // str listing selected applications
    )
    FUNFAM_HMMER_PARSER(FUNFAM_HMMER_RUNNER.out[0])  // hmmer.out pathS - one per cath gene3d superfam
    FUNFAM_CATH_RESOLVE_HITS(FUNFAM_HMMER_RUNNER.out)
    FUNFAM_FILTER_MATCHES(
        FUNFAM_HMMER_PARSER.out,           // add-superfams out file
        FUNFAM_CATH_RESOLVE_HITS.out,      // ips6 json
        FUNFAM_HMMER_RUNNER.out[1]         // post-processing-params
    )

    // HAMAP (+ pfsearch_wrapper.py)
    runner_hamap_params = fasta.combine(member_params.hamap)
    HAMAP_HMMER_RUNNER(runner_hamap_params)
    HAMAP_HMMER_PARSER(
        HAMAP_HMMER_RUNNER.out[0],  // hmmer.out path
    )
    HAMAP_POST_PROCESSER(
        HAMAP_HMMER_RUNNER.out[1],  // post-processing-params
        HAMAP_HMMER_RUNNER.out[2],  // path to fasta file
        HAMAP_HMMER_RUNNER.out[3],  // hmmer .tbl file path
    )
    HAMAP_FILTER_MATCHES(
        HAMAP_HMMER_PARSER.out,     // internal IPS6 JSON
        HAMAP_POST_PROCESSER.out,   // output from pfsearch_wrapper.py
    )

    // NCBIfam
    runner_hmmer_ncbifam_params = fasta.combine(member_params.ncbifam)
    NCBIFAM_HMMER_RUNNER(runner_hmmer_ncbifam_params)
    NCBIFAM_HMMER_PARSER(NCBIFAM_HMMER_RUNNER.out[0])  // hmmer.out path

    // Panther (+ treegrafter + epa-ng)
    runner_panther_params = fasta.combine(member_params.panther)
    PANTHER_HMMER_RUNNER(runner_panther_params)
    PANTHER_HMMER_PARSER(PANTHER_HMMER_RUNNER.out[0])  // hmmer.out path
    PANTHER_POST_PROCESSER(
        PANTHER_HMMER_RUNNER.out[0],  // hmmer.out path
        PANTHER_HMMER_RUNNER.out[1],  // post-processing-params
        fasta
    )
    PANTHER_FILTER_MATCHES(
        PANTHER_HMMER_PARSER.out,     // internal ips6 json
        PANTHER_POST_PROCESSER.out    // treegrafter output
    )

    // Pfam
    runner_hmmer_pfam_params = fasta.combine(member_params.pfam)
    PFAM_HMMER_RUNNER(runner_hmmer_pfam_params)
    PFAM_HMMER_PARSER(PFAM_HMMER_RUNNER.out[0])  // hmmer.out path
    PFAM_FILTER_MATCHES(
        PFAM_HMMER_PARSER.out,     // ips6 json
        PFAM_HMMER_RUNNER.out[1]   // post-processing-params
    )

    // PIRSF (+ filter_ips6_matches.py for post-processing)
    runner_pirsf_params = fasta.combine(member_params.pirsf)
    PIRSF_HMMER_RUNNER(runner_pirsf_params)
    PIRSF_HMMER_PARSER(PIRSF_HMMER_RUNNER.out[0])  // hmmer.out path
    PIRSF_FILTER_MATCHES(
        PIRSF_HMMER_PARSER.out,     // ips6 json
        PIRSF_HMMER_RUNNER.out[1],  // hmmer dtbl file -- needed to get tlen value
        PIRSF_HMMER_RUNNER.out[2]   // post-processing-params
    )

    // PIRSR
    runner_pirsr_params = fasta.combine(member_params.pirsr)
    PIRSR_HMMER_RUNNER(runner_pirsr_params)
    PIRSR_HMMER_PARSER(
        PIRSR_HMMER_RUNNER.out[0]  // out file
    )
    PIRSR_FILTER_MATCHES(
        PIRSR_HMMER_PARSER.out,  // ips6 json
        PIRSR_HMMER_RUNNER.out[1]  // post-processing-params
    )

    // SFLD (+ post-processing binary to add sites and filter hits)
    runner_sfld_params = fasta.combine(member_params.sfld)
    SFLD_HMMER_RUNNER(runner_sfld_params)
    SFLD_HMMER_PARSER(SFLD_HMMER_RUNNER.out[0])
    SFLD_POST_PROCESSER(SFLD_HMMER_RUNNER.out)   // hmmer.out, post-process params, alignment, dtbl file
    SFLD_FILTER_MATCHES(SFLD_HMMER_PARSER.out, SFLD_POST_PROCESSER.out)

    // SMART (HMMER2:hmmpfam + kinase filter)
    runner_smart_params = fasta.combine(member_params.smart)
    SMART_HMMER2_RUNNER(runner_smart_params)
    HMMER2_PARSER(SMART_HMMER2_RUNNER.out)
    SMART_FILTER_MATCHES(
        HMMER2_PARSER.out,
        SMART_HMMER2_RUNNER.out[1],
    )

    // Superfamily
    runner_hmmer_superfamily_params = fasta.combine(member_params.superfamily)
    SUPERFAMILY_HMMER_RUNNER(runner_hmmer_superfamily_params)
    SUPERFAMILY_POST_PROCESSER(
        SUPERFAMILY_HMMER_RUNNER.out[0],  // hmmer.out path
        SUPERFAMILY_HMMER_RUNNER.out[1],  // post-processing-params
        SUPERFAMILY_HMMER_RUNNER.out[2],  // fasta path
    )
    SUPERFAMILY_FILTER_MATCHES(
        SUPERFAMILY_POST_PROCESSER.out,
        SUPERFAMILY_HMMER_RUNNER.out[3],  // hmm path
    )

    /*
    Member databases that do NOT use HMMER
    */
    // CDD
    runner_cdd_params = fasta.combine(member_params.cdd)
    CDD_RUNNER(runner_cdd_params)
    CDD_POSTPROCESS(CDD_RUNNER.out)
    CDD_PARSER(CDD_POSTPROCESS.out)

    // COILS
    runner_coils_params = fasta.combine(member_params.coils)
    COILS_RUNNER(runner_coils_params)
    COILS_PARSER(COILS_RUNNER.out)

    // MOBIDB
    runner_mobidb_params = fasta.combine(member_params.mobidb)
    MOBIDB_RUNNER(runner_mobidb_params)
    MOBIDB_PARSER(MOBIDB_RUNNER.out)

    // PHOBIUS
    runner_phobius_params = fasta.combine(member_params.phobius)
    PHOBIUS_RUNNER(runner_phobius_params)
    PHOBIUS_PARSER(PHOBIUS_RUNNER.out)

    // PRINTS
    runner_prints_params = fasta.combine(member_params.prints)
    PRINTS_RUNNER(runner_prints_params)
    PRINTS_PARSER(PRINTS_RUNNER.out)

    // PROSITE Patterns (uses pfscanV3)
    runner_patterns = fasta.combine(member_params.prosite_patterns)
    PROSITE_PATTERNS_RUNNER(runner_patterns)
    PROSITE_PATTERNS_PARSER(PROSITE_PATTERNS_RUNNER.out)

    // PROSITE Profiles (uses pfsearchV3)
    runner_profiles = fasta.combine(member_params.prosite_profiles)
    PROSITE_PROFILES_RUNNER(runner_profiles)
    PROSITE_PROFILES_PARSER(PROSITE_PROFILES_RUNNER.out)

    // SignalP
    runner_signalp_params = fasta.combine(member_params.signalp)
    SIGNALP_RUNNER(runner_signalp_params)
    SIGNALP_PARSER(SIGNALP_RUNNER.out)

    // SignalP_euk
    runner_signalp_euk_params = fasta.combine(member_params.signalp_euk)
    SIGNALP_EUK_RUNNER(runner_signalp_euk_params)
    SIGNALP_EUK_PARSER(SIGNALP_EUK_RUNNER.out)

    /*
    Gather the results
    */
    if (applications.contains("gene3d")) {
        ANTIFAM_HMMER_PARSER.out.concat(
            NCBIFAM_HMMER_PARSER.out,
            GENE3D_FILTER_MATCHES.out[0],
            FUNFAM_FILTER_MATCHES.out,
            HAMAP_FILTER_MATCHES.out,
            PANTHER_FILTER_MATCHES.out,
            PFAM_FILTER_MATCHES.out,
            PIRSF_FILTER_MATCHES.out,
            PIRSR_FILTER_MATCHES.out,
            SFLD_FILTER_MATCHES.out,
            SMART_FILTER_MATCHES.out,
            CDD_PARSER.out,
            COILS_PARSER.out,
            MOBIDB_PARSER.out,
            PHOBIUS_PARSER.out,
            PRINTS_PARSER.out,
            PROSITE_PATTERNS_PARSER.out,
            PROSITE_PROFILES_PARSER.out,
            SIGNALP_PARSER.out,
            SIGNALP_EUK_PARSER.out,
            SUPERFAMILY_FILTER_MATCHES.out
        )
        .set { parsed_results }
    }
    else {
        ANTIFAM_HMMER_PARSER.out.concat(
            NCBIFAM_HMMER_PARSER.out,
            HAMAP_FILTER_MATCHES.out,
            FUNFAM_FILTER_MATCHES.out,
            PANTHER_FILTER_MATCHES.out,
            PFAM_FILTER_MATCHES.out,
            PIRSF_FILTER_MATCHES.out,
            PIRSR_FILTER_MATCHES.out,
            SFLD_FILTER_MATCHES.out,
            SMART_FILTER_MATCHES.out,
            CDD_PARSER.out,
            COILS_PARSER.out,
            MOBIDB_PARSER.out,
            PHOBIUS_PARSER.out,
            PRINTS_PARSER.out,
            PROSITE_PATTERNS_PARSER.out,
            PROSITE_PROFILES_PARSER.out,
            SIGNALP_PARSER.out,
            SIGNALP_EUK_PARSER.out,
            SUPERFAMILY_FILTER_MATCHES.out
        )
        .set { parsed_results }
    }

    emit:
    parsed_results
}
