include { ANTIFAM           } from "../antifam"
include { CATH              } from "../cath"
include { CDD               } from "../cdd"
include { COILS             } from "../coils"
include { DEEPTMHMM         } from "../deeptmhmm"
include { HAMAP             } from "../hamap"
include { INTERPRO_N        } from "../interpro-n"
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
include { TMBED             } from "../tmbed"
include { REPORT_NO_MATCHES } from "../../modules/no_matches"

workflow SCAN_SEQUENCES {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    db_releases         // map: [db: version, dirpath]
    applications        // list[str], applications to run in this workflow
    apps_config         // map of applications
    all_appls           // list[str], applications to run across all workflows 
    batch_size          // int, sub-batch size for computationally demanding apps

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        ANTIFAM(
            ch_seqs,
            db_releases.antifam.dirpath,
            apps_config.antifam.hmm
        )

        results = results.mix(ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        CATH(
            ch_seqs,
            applications.contains("cathgene3d"),
            db_releases.cathgene3d.dirpath,
            apps_config.cathgene3d.hmm,
            apps_config.cathgene3d.model2sfs,
            apps_config.cathgene3d.disc_regs,
            applications.contains("cathfunfam"),
            db_releases.cathfunfam.dirpath,
            batch_size
        ).set{ ch_cath }

        results = results.mix(ch_cath)
    }

    if (applications.contains("cdd")) {
        CDD(
            ch_seqs,
            db_releases.cdd.dirpath,
            apps_config.cdd.rpsblast_db,
            apps_config.cdd.rpsproc_db
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
            apps_config.deeptmhmm.dir,
            apps_config.deeptmhmm.use_gpu,
            batch_size
        )
        results = results.mix(DEEPTMHMM.out)
    }

    if (applications.contains("hamap")) {
        HAMAP(
            ch_seqs,
            db_releases.hamap.dirpath,
            apps_config.hamap.profiles
        )

        results = results.mix(HAMAP.out)
    }

    if (applications.contains("interpro_n")) {
        INTERPRO_N(
            ch_seqs,
            all_appls,
            apps_config.interpro_n.dir,
            apps_config.interpro_n.use_gpu,
            apps_config.interpro_n.batch_size,
            batch_size
        )

        results = results.mix(INTERPRO_N.out)
    }

    if (applications.contains("mobidblite")) {
        MOBIDBLITE(ch_seqs, batch_size)
        results = results.mix(MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        NCBIFAM(
            ch_seqs,
            db_releases.ncbifam.dirpath,
            apps_config.ncbifam.hmm
        )

        results = results.mix(NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        PANTHER(
            ch_seqs,
            db_releases.panther.dirpath,
            apps_config.panther.hmm,
            apps_config.panther.msf,
            batch_size
        )

        results = results.mix(PANTHER.out)
    }

    if (applications.contains("pfam")) {
        PFAM(
            ch_seqs,
            db_releases.pfam.dirpath,
            apps_config.pfam.hmm,
            apps_config.pfam.dat
        )

        results = results.mix(PFAM.out)
    }

    if (applications.contains("phobius")) {
        PHOBIUS(
            ch_seqs,
            apps_config.phobius.dir
        )
        results = results.mix(PHOBIUS.out)
    }

    if (applications.contains("pirsf")) {
        PIRSF(
            ch_seqs,
            db_releases.pirsf.dirpath,
            apps_config.pirsf.hmm,
            apps_config.pirsf.dat
        )

        results = results.mix(PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        PIRSR(
            ch_seqs,
            db_releases.pirsr.dirpath,
            apps_config.pirsr.hmm,
            apps_config.pirsr.rules
        )

        results = results.mix(PIRSR.out)
    }

    if (applications.contains("prints")) {
        PRINTS(
            ch_seqs,
            db_releases.prints.dirpath,
            apps_config.prints.pval,
            apps_config.prints.hierarchy,
            batch_size
        )

        results = results.mix(PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        PROSITE_PATTERNS(
            ch_seqs,
            db_releases.prositepatterns.dirpath,
            apps_config.prositepatterns.dat,
            apps_config.prositepatterns.evaluator
        )

        results = results.mix(PROSITE_PATTERNS.out)
    }

    if (applications.contains("prositeprofiles")) {
        PROSITE_PROFILES(
            ch_seqs,
            db_releases.prositeprofiles.dirpath,
            apps_config.prositeprofiles.profiles,
            apps_config.prositeprofiles.skip_flagged_profiles
        )

        results = results.mix(PROSITE_PROFILES.out)
    }

    if (applications.contains("sfld")) {
        SFLD(
            ch_seqs,
            db_releases.sfld.dirpath,
            apps_config.sfld.hmm,
            apps_config.sfld.sites_annotation,
            apps_config.sfld.hierarchy
        )

        results = results.mix(SFLD.out)
    }

    if (applications.contains("signalp_euk") || applications.contains("signalp_prok")) {
        SIGNALP(
            ch_seqs,
            applications,
            apps_config.signalp_euk.organism,
            apps_config.signalp_euk.mode,
            apps_config.signalp_euk.dir,
            apps_config.signalp_euk.use_gpu,
            apps_config.signalp_prok.organism,
            apps_config.signalp_prok.mode,
            apps_config.signalp_prok.dir,
            apps_config.signalp_prok.use_gpu,
            batch_size
        ).set{ ch_signalp }
        results = results.mix(ch_signalp)
    }

    if (applications.contains("smart")) {
        SMART(
            ch_seqs,
            db_releases.smart.dirpath,
            apps_config.smart.hmm3,
            apps_config.smart.hmm2,
            apps_config.smart.chunk_size
        )

        results = results.mix(SMART.out)
    }

    if (applications.contains("superfamily")) {
        SUPERFAMILY(
            ch_seqs,
            db_releases.superfamily.dirpath,
            apps_config.superfamily.hmm,
            apps_config.superfamily.selfhits,
            apps_config.superfamily.cla,
            apps_config.superfamily.model,
            apps_config.superfamily.pdbj95d,
            batch_size
        )

        results = results.mix(SUPERFAMILY.out)
    }

    if (applications.contains("tmbed")) {
        TMBED(
            ch_seqs,
            apps_config.tmbed.use_gpu,
            apps_config.tmbed.chunk_size,
            apps_config.tmbed.chunk_overlap,
            apps_config.tmbed.smooth_window,
            apps_config.tmbed.batch_size,
            batch_size
        )

        results = results.mix(TMBED.out)
    }

    ch_results = ch_seqs.join(
        results.groupTuple(),
        failOnDuplicate: false, // may happen when applications perfom sub-batching
        failOnMismatch: false   // may happen when all sequences found in Matches API
    )

    ch_no_matches = REPORT_NO_MATCHES(ch_results)

    merged_results = ch_results
        .join(
            ch_no_matches,
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { meta, fasta, member_paths, no_match_path ->
            [meta, member_paths + [no_match_path]]
        }

    emit:
    merged_results
}