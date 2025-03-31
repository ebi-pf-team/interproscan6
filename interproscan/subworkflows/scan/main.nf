include { RUN_ANTIFAM; PARSE_ANTIFAM                                              } from  "../../modules/antifam"
include { SEARCH_GENE3D; RESOLVE_GENE3D; ASSIGN_CATH; PARSE_CATHGENE3D            } from  "../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; RESOLVE_FUNFAM; PARSE_FUNFAM             } from  "../../modules/cath/funfam"
include { RUN_RPSBLAST; RUN_RPSPROC; PARSE_RPSPROC                                } from  "../../modules/cdd"
include { RUN_COILS; PARSE_COILS                                                  } from  "../../modules/coils"
include { PREPROCESS_HAMAP; PREPARE_HAMAP; RUN_HAMAP; PARSE_HAMAP                 } from  "../../modules/hamap"
include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE                                        } from  "../../modules/mobidblite"
include { RUN_NCBIFAM; PARSE_NCBIFAM                                              } from  "../../modules/ncbifam"
include { SEARCH_PANTHER; PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER     } from  "../../modules/panther"
include { SEARCH_PFAM; PARSE_PFAM                                                 } from  "../../modules/pfam"
include { SEARCH_PHOBIUS; PARSE_PHOBIUS                                           } from  "../../modules/phobius"
include { RUN_PIRSR; PARSE_PIRSR                                                  } from  "../../modules/pirsr"
include { RUN_PIRSF; PARSE_PIRSF                                                  } from  "../../modules/pirsf"
include { RUN_PRINTS; PARSE_PRINTS                                                } from  "../../modules/prints"
include { RUN_PFSCAN ; PARSE_PFSCAN                                               } from  "../../modules/prosite/patterns"
include { RUN_PFSEARCH ; PARSE_PFSEARCH                                           } from  "../../modules/prosite/profiles"
include { RUN_SFLD; POST_PROCESS_SFLD; PARSE_SFLD                                 } from  "../../modules/sfld"
include { RUN_SIGNALP as RUN_SIGNALP_EUK; PARSE_SIGNALP as PARSE_SIGNALP_EUK      } from  "../../modules/signalp"
include { RUN_SIGNALP as RUN_SIGNALP_PROK; PARSE_SIGNALP as PARSE_SIGNALP_PROK    } from  "../../modules/signalp"
include { SEARCH_SMART; PARSE_SMART                                               } from  "../../modules/smart"
include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY                                   } from  "../../modules/superfamily"
include { RUN_DEEPTMHMM; PARSE_DEEPTMHMM                                          } from  "../../modules/tmhmm"

workflow SCAN_SEQUENCES {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    member_db_releases  // map: [db: version number]
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        def antifam_release = member_db_releases['antifam']
        def antifam_dir = "${datadir}/${appsConfig.antifam.dir}/${antifam_release}"
        RUN_ANTIFAM(
            ch_seqs,
            "${antifam_dir}/${appsConfig.antifam.hmm}"
        )

        PARSE_ANTIFAM(RUN_ANTIFAM.out)
        results = results.mix(PARSE_ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        def gene3d_release = member_db_releases['cath-gene3d']
        def gene3d_dir = "${datadir}/${appsConfig.cathgene3d.dir}/${gene3d_release}"
        // Search Gene3D profiles
        SEARCH_GENE3D(
            ch_seqs,
            "${gene3d_dir}/${appsConfig.cathgene3d.hmm}")
        ch_gene3d = SEARCH_GENE3D.out

        // Select best domain matches
        RESOLVE_GENE3D(ch_gene3d)

        // Assign CATH superfamily to matches
        ASSIGN_CATH(
            RESOLVE_GENE3D.out,
            "${gene3d_dir}/${appsConfig.cathgene3d.model2sfs}",
            "${gene3d_dir}/${appsConfig.cathgene3d.disc_regs}"
        )

        // Join results and parse them
        ch_cathgene3d = ch_gene3d.join(ASSIGN_CATH.out)
        PARSE_CATHGENE3D(ch_cathgene3d)

        if (applications.contains("cathgene3d")) {
            results = results.mix(PARSE_CATHGENE3D.out)
        }

        if (applications.contains("cathfunfam")) {
            def funfam_release = member_db_releases['cath-funfam']
            def funfam_dir = "${datadir}/${appsConfig.cathfunfam.dir}/${funfam_release}"
            // Find unique CATH superfamilies with at least one hit
            PREPARE_FUNFAM(
                PARSE_CATHGENE3D.out,
                "${funfam_dir}"
            )

            /*
                Join input fasta file with superfamilies.
                We split in smaller chunks to parallelize searching FunFam profiles
            */
            ch_seqs
                .join(PREPARE_FUNFAM.out)
                .flatMap { id, file, supfams ->
                    supfams
                        .collate(appsConfig.cathfunfam.chunkSize)
                        .collect { chunk -> [id, file, chunk] }
                }
                .set { ch_funfams }

            // Search FunFam profiles
            SEARCH_FUNFAM(
                ch_funfams,
                "${funfam_dir}"
            )

            // Select best domain matches
            RESOLVE_FUNFAM(SEARCH_FUNFAM.out)

            // Join results and parse them
            ch_cathfunfam = SEARCH_FUNFAM.out.join(RESOLVE_FUNFAM.out)
            PARSE_FUNFAM(ch_cathfunfam)
            results = results.mix(PARSE_FUNFAM.out)
        }
    }

    if (applications.contains("cdd")) {
        def cdd_release = member_db_releases['cdd']
        def cdd_dir = "${datadir}/${appsConfig.cdd.dir}/${cdd_release}"
        RUN_RPSBLAST(
            ch_seqs,
            "${cdd_dir}/${appsConfig.cdd.rpsblast_db}"
        )

        RUN_RPSPROC(
            RUN_RPSBLAST.out,
            "${cdd_dir}/${appsConfig.cdd.rpsproc_db}"
        )

        PARSE_RPSPROC(RUN_RPSPROC.out)
        results = results.mix(PARSE_RPSPROC.out)
    }

    if (applications.contains("coils")) {
        RUN_COILS(ch_seqs)
        PARSE_COILS(RUN_COILS.out)
        results = results.mix(PARSE_COILS.out)
    }

    if (applications.contains("deeptmhmm")) {
        RUN_DEEPTMHMM(
            ch_seqs,
            appsConfig.deeptmhmm.dir
        )
        PARSE_DEEPTMHMM(RUN_DEEPTMHMM.out)
        results = results.mix(PARSE_DEEPTMHMM.out)
    }

    if (applications.contains("hamap")) {
        def hamap_release = member_db_releases['hamap']
        def hamap_dir = "${datadir}/${appsConfig.hamap.dir}/${hamap_release}"
        PREPROCESS_HAMAP(
            ch_seqs,
            "${hamap_dir}/${appsConfig.hamap.hmm}"
        )

        PREPARE_HAMAP(
            PREPROCESS_HAMAP.out,
            "${hamap_dir}/${appsConfig.hamap.dir}"
        )

        RUN_HAMAP(
            PREPARE_HAMAP.out,
            "${hamap_dir}/${appsConfig.hamap.dir}"
        )

        PARSE_HAMAP(RUN_HAMAP.out)
        results = results.mix(PARSE_HAMAP.out)
    }

    if (applications.contains("mobidblite")) {
        RUN_MOBIDBLITE(ch_seqs)

        PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out)
        results = results.mix(PARSE_MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        def ncbifam_release = member_db_releases['ncbifam']
        def ncbifam_dir = "${datadir}/${appsConfig.ncbifam.dir}/${ncbifam_release}"
        RUN_NCBIFAM(
            ch_seqs,
            "${ncbifam_dir}/${appsConfig.ncbifam.hmm}"
        )

        PARSE_NCBIFAM(RUN_NCBIFAM.out)
        results = results.mix(PARSE_NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        def panther_release = member_db_releases['panther']
        def panther_dir = "${datadir}/${appsConfig.panther.dir}/${panther_release}"
        SEARCH_PANTHER(
            ch_seqs,
            "${panther_dir}/${appsConfig.panther.hmm}"
        )
        ch_panther = SEARCH_PANTHER.out

        PREPARE_TREEGRAFTER(
            ch_panther,
            "${panther_dir}/${appsConfig.panther.msf}"
        )

        RUN_TREEGRAFTER(
            PREPARE_TREEGRAFTER.out.fasta,
            "${panther_dir}/${appsConfig.panther.msf}"
        )

        PARSE_PANTHER(
            PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out)
        )
        results = results.mix(PARSE_PANTHER.out)
    }

    if (applications.contains("phobius")) {
        SEARCH_PHOBIUS(
            ch_seqs,
            appsConfig.phobius.dir
        )

        PARSE_PHOBIUS(SEARCH_PHOBIUS.out)
        results = results.mix(PARSE_PHOBIUS.out)
    }

    if (applications.contains("pfam")) {
        def pfam_release = member_db_releases['pfam']
        def pfam_dir = "${datadir}/${appsConfig.pfam.dir}/${pfam_release}"
        SEARCH_PFAM(
            ch_seqs,
            "${pfam_dir}/${appsConfig.pfam.hmm}"
        )

        PARSE_PFAM(
            SEARCH_PFAM.out,
            "${pfam_dir}/${appsConfig.pfam.seed}",
            "${pfam_dir}/${appsConfig.pfam.clan}",
            "${pfam_dir}/${appsConfig.pfam.dat}"
        )
        results = results.mix(PARSE_PFAM.out)
    }

    if (applications.contains("pirsf")) {
        def pirsf_release = member_db_releases['pirsf']
        def pirsf_dir = "${datadir}/${appsConfig.pirsf.dir}/${pirsf_release}"
        RUN_PIRSF(ch_seqs,
            "${pirsf_dir}/${appsConfig.pirsf.hmm}"
        )
        PARSE_PIRSF(RUN_PIRSF.out,
            "${pirsf_dir}/${appsConfig.pirsf.dat}")

        results = results.mix(PARSE_PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        def pirsr_release = member_db_releases['pirsr']
        def pirsr_dir = "${datadir}/${appsConfig.pirsr.dir}/${pirsr_release}"
        RUN_PIRSR(ch_seqs,
            "${pirsr_dir}/${appsConfig.pirsr.hmm}")

        PARSE_PIRSR(RUN_PIRSR.out,
            "${pirsr_dir}/${appsConfig.pirsr.rules}")

        results = results.mix(PARSE_PIRSR.out)
    }

    if (applications.contains("prints")) {
        def prints_release = member_db_releases['prints']
        def prints_dir = "${datadir}/${appsConfig.prints.dir}/${prints_release}"
        RUN_PRINTS(
            ch_seqs,
            "${prints_dir}/${appsConfig.prints.pval}"
        )

        PARSE_PRINTS(
            RUN_PRINTS.out,
            "${prints_dir}/${appsConfig.prints.hierarchy}"
        )
        results = results.mix(PARSE_PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        def prositepatterns_release = member_db_releases['prosite patterns']
        def prositepatterns_dir = "${datadir}/${appsConfig.prositepatterns.dir}/${prositepatterns_release}"
        RUN_PFSCAN(
            ch_seqs,
            "${prositepatterns_dir}/${appsConfig.prositepatterns.dat}",
            "${prositepatterns_dir}/${appsConfig.prositepatterns.evaluator}")

        PARSE_PFSCAN(RUN_PFSCAN.out)
        results = results.mix(PARSE_PFSCAN.out)
    }

    if (applications.contains("prositeprofiles")) {
        def prositeprofiles_release = member_db_releases['prosite profiles']
        def prositeprofiles_dir = "${datadir}/${appsConfig.prositeprofiles.dir}/${prositeprofiles_release}"
        RUN_PFSEARCH(
            ch_seqs,
            "${prositeprofiles_dir}/${appsConfig.prositeprofiles.dir}")
        PARSE_PFSEARCH(
            RUN_PFSEARCH.out,
            "${prositeprofiles_dir}/${appsConfig.prositeprofiles.skip_flagged_profiles}")
        results = results.mix(PARSE_PFSEARCH.out)
    }

    if (applications.contains("sfld")) {
        def sfld_release = member_db_releases['sfld']
        def sfld_dir = "${datadir}/${appsConfig.sfld.dir}/${sfld_release}"
        RUN_SFLD(ch_seqs,
            "${sfld_dir}/${appsConfig.sfld.hmm}")

        POST_PROCESS_SFLD(RUN_SFLD.out,
            "${sfld_dir}/${appsConfig.sfld.sites_annotation}")

        PARSE_SFLD(POST_PROCESS_SFLD.out,
            "${sfld_dir}/${appsConfig.sfld.hierarchy}")

        results = results.mix(PARSE_SFLD.out)
    }

    if (applications.contains("signalp_euk")) {
        RUN_SIGNALP_EUK(
            ch_seqs,
            appsConfig.signalp_euk.organism,
            appsConfig.signalp_euk.mode,
            appsConfig.signalp_euk.dir
        )
        PARSE_SIGNALP_EUK(RUN_SIGNALP_EUK.out)
        results = results.mix(PARSE_SIGNALP_EUK.out)
    }

    if (applications.contains("signalp_prok")) {
        RUN_SIGNALP_PROK(
            ch_seqs,
            appsConfig.signalp_prok.organism,
            appsConfig.signalp_prok.mode,
            appsConfig.signalp_prok.dir
        )
        PARSE_SIGNALP_PROK(RUN_SIGNALP_PROK.out)
        results = results.mix(PARSE_SIGNALP_PROK.out)
    }

    if (applications.contains("smart")) {
        def smart_release = member_db_releases['smart']
        def smart_dir = "${datadir}/${appsConfig.smart.dir}/${smart_release}"
        SEARCH_SMART(
            ch_seqs,
            "${smart_dir}/${appsConfig.smart.hmmbin}"
        )

        PARSE_SMART(
            SEARCH_SMART.out,
            "${smart_dir}/${appsConfig.smart.hmm}"
        )
        results = results.mix(PARSE_SMART.out)
    }

    if (applications.contains("superfamily")) {
        def superfamily_release = member_db_releases['superfamily']
        def superfamily_dir = "${datadir}/${appsConfig.superfamily.dir}/${superfamily_release}"
        SEARCH_SUPERFAMILY(
            ch_seqs,
            "${superfamily_dir}/${appsConfig.superfamily.hmm}",
            "${superfamily_dir}/${appsConfig.superfamily.selfhits}",
            "${superfamily_dir}/${appsConfig.superfamily.cla}",
            "${superfamily_dir}/${appsConfig.superfamily.model}",
            "${superfamily_dir}/${appsConfig.superfamily.pdbj95d}"
        )

        PARSE_SUPERFAMILY(
            SEARCH_SUPERFAMILY.out,
            "${superfamily_dir}/${appsConfig.superfamily.model}",
            "${superfamily_dir}/${appsConfig.superfamily.hmm}"
        )
        results = results.mix(PARSE_SUPERFAMILY.out)
    }

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}