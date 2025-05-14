include { ANTIFAM           } from "../antifam"
include { CATH              } from "../cath"
include { CDD               } from "../cdd"
include { COILS             } from "../coils"
include { DEEPTMHMM         } from "../deeptmhmm"
include { HAMAP             } from "../hamap"
include { MOBIDBLITE        } from "../mobidblite"
include { NCBIFAM           } from "../ncbifam"
include { PANTHER           } from "../panther"
include { PFAM              } from "../pfam"
include { PHOBIUS           } from "../phobius"
include { PIRSF             } from "../pirsf"
include { PIRSR             } from "../pirsr"
include { PRINTS            } from "../prints"
include { PROSITE_PATTERNS  } from "../prosite/patterns"
include { PROSITE_PROFILES  } from "../prosite/profiles"
include { SFLD              } from "../sfld"
include { SIGNALP           } from "../signalp"
include { SMART             } from "../smart"
include { SUPERFAMILY       } from "../superfamily"
include { REPORT_NO_MATCHES } from "../../modules/no_matches"

workflow SCAN_SEQUENCES {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    db_releases         // map: [db: version, dirpath]
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory
    signalp_gpu
    deeptmhmm_gpu

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        ANTIFAM(
            ch_seqs,
            db_releases.antifam.dirpath,
            appsConfig.antifam.hmm
        )

        results = results.mix(ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        CATH(
            ch_seqs,
            applications.contains("cathgene3d"),
            db_releases.cathgene3d.dirpath,
            appsConfig.cathgene3d.hmm,
            appsConfig.cathgene3d.model2sfs,
            appsConfig.cathgene3d.disc_regs,
            applications.contains("cathfunfam"),
            db_releases.cathfunfam.dirpath
        ).set{ ch_cath }

        results = results.mix(ch_cath)
    }

    if (applications.contains("cdd")) {
        CDD(
            ch_seqs,
            db_releases.cdd.dirpath,
            appsConfig.cdd.rpsblast_db,
            appsConfig.cdd.rpsproc_db
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
            db_releases.hamap.dirpath,
            appsConfig.hamap.hmm,
            appsConfig.hamap.profiles
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
            db_releases.ncbifam.dirpath,
            appsConfig.ncbifam.hmm
        )

        results = results.mix(NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        PANTHER(
            ch_seqs,
            db_releases.panther.dirpath,
            appsConfig.panther.hmm,
            appsConfig.panther.msf
        )

        results = results.mix(PANTHER.out)
    }

    if (applications.contains("pfam")) {
        PFAM(
            ch_seqs,
            db_releases.pfam.dirpath,
            appsConfig.pfam.hmm,
            appsConfig.pfam.dat
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
            db_releases.pirsf.dirpath,
            appsConfig.pirsf.hmm,
            appsConfig.pirsf.dat
        )

        results = results.mix(PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        PIRSR(
            ch_seqs,
            db_releases.pirsr.dirpath,
            appsConfig.pirsr.hmm,
            appsConfig.pirsr.rules
        )

        results = results.mix(PIRSR.out)
    }

    if (applications.contains("prints")) {
        PRINTS(
            ch_seqs,
            db_releases.prints.dirpath,
            appsConfig.prints.pval,
            appsConfig.prints.hierarchy
        )

        results = results.mix(PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        PROSITE_PATTERNS(
            ch_seqs,
            db_releases.prositepatterns.dirpath,
            appsConfig.prositepatterns.dat,
            appsConfig.prositepatterns.evaluator
        )

        results = results.mix(PROSITE_PATTERNS.out)
    }

    if (applications.contains("prositeprofiles")) {
        PROSITE_PROFILES(
            ch_seqs,
            db_releases.prositeprofiles.dirpath,
            appsConfig.prositeprofiles.profiles,
            appsConfig.prositeprofiles.skip_flagged_profiles
        )

        results = results.mix(PROSITE_PROFILES.out)
    }

    if (applications.contains("sfld")) {
        SFLD(
            ch_seqs,
            db_releases.sfld.dirpath,
            appsConfig.sfld.hmm,
            appsConfig.sfld.sites_annotation,
            appsConfig.sfld.hierarchy
        )

        results = results.mix(SFLD.out)
    }

    if (applications.contains("signalp_euk") || applications.contains("signalp_prok")) {
        SIGNALP(
            ch_seqs,
            applications,
            appsConfig.signalp_euk.organism,
            appsConfig.signalp_euk.mode,
            appsConfig.signalp_euk.dir,
            appsConfig.signalp_prok.organism,
            appsConfig.signalp_prok.mode,
            appsConfig.signalp_prok.dir,
            signalp_gpu
        ).set{ ch_signalp }
        results = results.mix(ch_signalp)
    }

    if (applications.contains("smart")) {
        SMART(
            ch_seqs,
            db_releases.smart.dirpath,
            appsConfig.smart.hmm3,
            appsConfig.smart.hmm2,
            appsConfig.smart.chunk_size
        )

        results = results.mix(SMART.out)
    }

    if (applications.contains("superfamily")) {
        SUPERFAMILY(
            ch_seqs,
            db_releases.superfamily.dirpath,
            appsConfig.superfamily.hmm,
            appsConfig.superfamily.selfhits,
            appsConfig.superfamily.cla,
            appsConfig.superfamily.model,
            appsConfig.superfamily.pdbj95d
        )

        results = results.mix(SUPERFAMILY.out)
    }

    results
        .groupTuple()
        .set { grouped_results }

    ch_no_matches = REPORT_NO_MATCHES(grouped_results, ch_seqs)

    merged_results = grouped_results
        .join(ch_no_matches)
        .map { batch_idx, paths_list, no_match_path ->
            [batch_idx, paths_list + [no_match_path]]
        }

    emit:
    merged_results
}