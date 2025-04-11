include { RUN_PFSEARCH ; PARSE_PFSEARCH                                           } from  "../../modules/prosite/profiles"
include { RUN_SFLD; POST_PROCESS_SFLD; PARSE_SFLD                                 } from  "../../modules/sfld"
include { RUN_SIGNALP as RUN_SIGNALP_EUK; PARSE_SIGNALP as PARSE_SIGNALP_EUK      } from  "../../modules/signalp"
include { RUN_SIGNALP as RUN_SIGNALP_PROK; PARSE_SIGNALP as PARSE_SIGNALP_PROK    } from  "../../modules/signalp"
include { SCAN_SMART; PREPARE_SMART; SEARCH_SMART; PARSE_SMART                    } from  "../../modules/smart"
include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY                                   } from  "../../modules/superfamily"

include { ANTIFAM          } from "../applications/antifam"
include { CATH             } from "../applications/cath"
include { CDD              } from "../applications/cdd"
include { COILS            } from "../applications/coils"
include { DEEPTMHMM        } from "../applications/deeptmhmm"
include { HAMAP            } from "../applications/hamap"
include { MOBIDBLITE       } from "../applications/mobidblite"
include { NCBIFAM          } from "../applications/ncbifam"
include { PANTHER          } from "../applications/panther"
include { PFAM             } from "../applications/pfam"
include { PHOBIUS          } from "../applications/phobius"
include { PIRSF            } from "../applications/pirsf"
include { PIRSR            } from "../applications/pirsr"
include { PRINTS           } from "../applications/prints"
include { PROSITE_PATTERNS } from "../applications/prosite/patterns"
// include { PROSITE_PROFILE  } from "../applications/prosite/profiles"

workflow SCAN_SEQUENCES {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        ANTIFAM(
            ch_seqs,
            "${datadir}/${appsConfig.antifam.hmm}"
        )
        results = results.mix(ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        CATH(
            ch_seqs,
            applications,
            "${datadir}/${appsConfig.cathgene3d.hmm}",
            "${datadir}/${appsConfig.cathgene3d.model2sfs}",
            "${datadir}/${appsConfig.cathgene3d.disc_regs}",
            "${datadir}/${appsConfig.cathfunfam.dir}",
            appsConfig.cathfunfam.chunkSize
        ).set{ ch_cath }
        results = results.mix(ch_cath)
    }

    if (applications.contains("cdd")) {
        CDD(
            ch_seqs,
            "${datadir}/${appsConfig.cdd.rpsblast_db}",
            "${datadir}/${appsConfig.cdd.rpsproc_db}"
        )
        results = results.mix(CDD.out)
    }

    if (applications.contains("coils")) {
        COILS(ch_seqs)
        results = results.mix(COILS.out)
    }

    if (applications.contains("deeptmhmm")) {
        DEEPTMHMM(
            ch_seqs,
            appsConfig.deeptmhmm.dir
        )
        results = results.mix(DEEPTMHMM.out)
    }

    if (applications.contains("hamap")) {
        HAMAP(
            ch_seqs,
            "${datadir}/${appsConfig.hamap.hmm}",
            "${datadir}/${appsConfig.hamap.dir}"
        )
        results = results.mix(HAMAP.out)
    }

    if (applications.contains("mobidblite")) {
        MOBIDBLITE(ch_seqs)
        results = results.mix(MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        NCBIFAM(
            ch_seqs,
            "${datadir}/${appsConfig.ncbifam.hmm}"
        )
        results = results.mix(NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        PANTHER(
            ch_seqs,
            "${datadir}/${appsConfig.panther.hmm}",
            "${datadir}/${appsConfig.panther.msf}"
        )
        results = results.mix(PANTHER.out)
    }

    if (applications.contains("pfam")) {
        PFAM(
            ch_seqs,
            "${datadir}/${appsConfig.pfam.hmm}",
            "${datadir}/${appsConfig.pfam.dat}"
        )
        results = results.mix(PFAM.out)
    }

    if (applications.contains("phobius")) {
        PHOBIUS(
            ch_seqs,
            appsConfig.phobius.dir
        )
        results = results.mix(PHOBIUS.out)
    }

    if (applications.contains("pirsf")) {
        PIRSF(
            ch_seqs,
            "${datadir}/${appsConfig.pirsf.hmm}",
            "${datadir}/${appsConfig.pirsf.dat}"
        )
        results = results.mix(PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        PIRSR(
            ch_seqs,
            "${datadir}/${appsConfig.pirsr.hmm}",
            "${datadir}/${appsConfig.pirsr.rules}"
        )
        results = results.mix(PIRSR.out)
    }

    if (applications.contains("prints")) {
        PRINTS(
            ch_seqs,
            "${datadir}/${appsConfig.prints.pval}",
            "${datadir}/${appsConfig.prints.hierarchy}"
        )
        results = results.mix(PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        PROSITE_PATTERNS(
            ch_seqs,
            "${datadir}/${appsConfig.prositepatterns.dat}",
            "${datadir}/${appsConfig.prositepatterns.evaluator}"
        )
        results = results.mix(PROSITE_PATTERNS.out)
    }

    if (applications.contains("prositeprofiles")) {
        PROSITE_PROFILES(
            ch_seqs,
            "${datadir}/${appsConfig.prositeprofiles.dir}",
            "${datadir}/${appsConfig.prositeprofiles.skip_flagged_profiles}"
        )
        results = results.mix(PARSE_PFSEARCH.out)
    }

    if (applications.contains("sfld")) {
        RUN_SFLD(ch_seqs,
            "${datadir}/${appsConfig.sfld.hmm}")

        POST_PROCESS_SFLD(RUN_SFLD.out,
            "${datadir}/${appsConfig.sfld.sites_annotation}")

        PARSE_SFLD(POST_PROCESS_SFLD.out,
            "${datadir}/${appsConfig.sfld.hierarchy}")

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
        SCAN_SMART(
            ch_seqs,
            "${datadir}/${appsConfig.smart.hmm}"
        )

        PREPARE_SMART(
            SCAN_SMART.out,
            appsConfig.smart.chunkSize,
            "${datadir}/${appsConfig.smart.hmm_dir}"
        )

        SEARCH_SMART(
            PREPARE_SMART.out,
            "${datadir}/${appsConfig.smart.hmm_dir}"
        )

        PARSE_SMART(
            SEARCH_SMART.out,
            "${datadir}/${appsConfig.smart.hmm}"
        )
        results = results.mix(PARSE_SMART.out)
    }

    if (applications.contains("superfamily")) {
        SEARCH_SUPERFAMILY(
            ch_seqs,
            "${datadir}/${appsConfig.superfamily.hmm}",
            "${datadir}/${appsConfig.superfamily.selfhits}",
            "${datadir}/${appsConfig.superfamily.cla}",
            "${datadir}/${appsConfig.superfamily.model}",
            "${datadir}/${appsConfig.superfamily.pdbj95d}"
        )

        PARSE_SUPERFAMILY(
            SEARCH_SUPERFAMILY.out,
            "${datadir}/${appsConfig.superfamily.model}",
            "${datadir}/${appsConfig.superfamily.hmm}"
        )
        results = results.mix(PARSE_SUPERFAMILY.out)
    }

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}