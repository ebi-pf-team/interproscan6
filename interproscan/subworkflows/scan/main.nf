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
include { PROSITE_PROFILES } from "../applications/prosite/profiles"
include { SFLD             } from "../applications/sfld"
include { SIGNALP          } from "../applications/signalp"
include { SMART            } from "../applications/smart"
include { SUPERFAMILY      } from "../applications/superfamily"

workflow SCAN_SEQUENCES {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    databases           // map: [db: version, dirpath]
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        def antifam_dir = databases["antifam"]["dirpath"]

        ANTIFAM(
            ch_seqs,
            "${antifam_dir}/${appsConfig.antifam.hmm}"
        )

        results = results.mix(ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        def gene3d_release = member_db_releases['cath-gene3d']
        def gene3d_dir = "${datadir}/${appsConfig.cathgene3d.dir}/${gene3d_release}"

        def funfam_release = member_db_releases['cath-funfam']
        def funfam_dir = "${datadir}/${appsConfig.cathfunfam.dir}/${funfam_release}"

        CATH(
            ch_seqs,
            applications,
            "${gene3d_dir}/${appsConfig.cathgene3d.hmm}",
            "${gene3d_dir}/${appsConfig.cathgene3d.model2sfs}",
            "${gene3d_dir}/${appsConfig.cathgene3d.disc_regs}",
            "${funfam_dir}/${appsConfig.cathfunfam.dir}",
            appsConfig.cathfunfam.chunkSize
        ).set{ ch_cath }

        results = results.mix(ch_cath)
    }

    if (applications.contains("cdd")) {
        def cdd_release = member_db_releases['cdd']
        def cdd_dir = "${datadir}/${appsConfig.cdd.dir}/${cdd_release}"

        CDD(
            ch_seqs,
            "${cdd_dir}/${appsConfig.cdd.rpsblast_db}",
            "${cdd_dir}/${appsConfig.cdd.rpsproc_db}"
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
        def hamap_release = member_db_releases['hamap']
        def hamap_dir = "${datadir}/${appsConfig.hamap.dir}/${hamap_release}"

        HAMAP(
            ch_seqs,
            "${hamap_dir}/${appsConfig.hamap.hmm}",
            "${hamap_dir}/${appsConfig.hamap.dir}"
        )

        results = results.mix(HAMAP.out)
    }

    if (applications.contains("mobidblite")) {
        MOBIDBLITE(ch_seqs)
        results = results.mix(MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        def ncbifam_release = member_db_releases['ncbifam']
        def ncbifam_dir = "${datadir}/${appsConfig.ncbifam.dir}/${ncbifam_release}"
        NCBIFAM(
            ch_seqs,
            "${ncbifam_dir}/${appsConfig.ncbifam.hmm}"
        )

        results = results.mix(NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        def panther_release = member_db_releases['panther']
        def panther_dir = "${datadir}/${appsConfig.panther.dir}/${panther_release}"

        PANTHER(
            ch_seqs,
            "${panther_dir}/${appsConfig.panther.hmm}",
            "${panther_dir}/${appsConfig.panther.msf}"
        )

        results = results.mix(PANTHER.out)
    }

    if (applications.contains("pfam")) {
        def pfam_release = member_db_releases['pfam']
        def pfam_dir = "${datadir}/${appsConfig.pfam.dir}/${pfam_release}"

        PFAM(
            ch_seqs,
            "${pfam_dir}/${appsConfig.pfam.hmm}",
            "${pfam_dir}/${appsConfig.pfam.dat}"
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
        def pirsf_release = member_db_releases['pirsf']
        def pirsf_dir = "${datadir}/${appsConfig.pirsf.dir}/${pirsf_release}"
    
        PIRSF(
            ch_seqs,
            "${pirsf_dir}/${appsConfig.pirsf.hmm}",
            "${pirsf_dir}/${appsConfig.pirsf.dat}"
        )

        results = results.mix(PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        def pirsr_release = member_db_releases['pirsr']
        def pirsr_dir = "${datadir}/${appsConfig.pirsr.dir}/${pirsr_release}"

        PIRSR(
            ch_seqs,
            "${pirsr_dir}/${appsConfig.pirsr.hmm}",
            "${pirsr_dir}/${appsConfig.pirsr.rules}"
        )

        results = results.mix(PIRSR.out)
    }

    if (applications.contains("prints")) {
        def prints_release = member_db_releases['prints']
        def prints_dir = "${datadir}/${appsConfig.prints.dir}/${prints_release}"

        PRINTS(
            ch_seqs,
            "${prints_dir}/${appsConfig.prints.pval}",
            "${prints_dir}/${appsConfig.prints.hierarchy}"
        )

        results = results.mix(PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        def prositepatterns_release = member_db_releases['prosite patterns']
        def prositepatterns_dir = "${datadir}/${appsConfig.prositepatterns.dir}/${prositepatterns_release}"

        PROSITE_PATTERNS(
            ch_seqs,
            "${prositepatterns_dir}/${appsConfig.prositepatterns.dat}",
            "${prositepatterns_dir}/${appsConfig.prositepatterns.evaluator}"
        )

        results = results.mix(PROSITE_PATTERNS.out)
    }

    if (applications.contains("prositeprofiles")) {
        def prositeprofiles_release = member_db_releases['prosite profiles']
        def prositeprofiles_dir = "${datadir}/${appsConfig.prositeprofiles.dir}/${prositeprofiles_release}"

        PROSITE_PROFILES(
            ch_seqs,
            "${prositeprofiles_dir}/${appsConfig.prositeprofiles.dir}",
            "${prositeprofiles_dir}/${appsConfig.prositeprofiles.skip_flagged_profiles}"
        )

        results = results.mix(PROSITE_PROFILES.out)
    }

    if (applications.contains("sfld")) {
        def sfld_release = member_db_releases['sfld']
        def sfld_dir = "${datadir}/${appsConfig.sfld.dir}/${sfld_release}"

        SFLD(
            ch_seqs,
            "${sfld_dir}/${appsConfig.sfld.hmm}",
            "${sfld_dir}/${appsConfig.sfld.sites_annotation}",
            "${sfld_dir}/${appsConfig.sfld.hierarchy}"
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
            appsConfig.signalp_prok.dir
        ).set{ ch_signalp }
        results = results.mix(ch_signalp)
    }

    if (applications.contains("smart")) {
        def smart_release = member_db_releases['smart']
        def smart_dir = "${datadir}/${appsConfig.smart.dir}/${smart_release}"

        SMART(
            ch_seqs,
            "${smart_dir}/${appsConfig.smart.hmmer3_hmm}",
            "${smart_dir}/${appsConfig.smart.hmmer2_hmm}",
            "${smart_dir}/${appsConfig.smart.hmm_dir}",
            appsConfig.smart.chunkSize
        )

        results = results.mix(SMART.out)
    }

    if (applications.contains("superfamily")) {
        def superfamily_release = member_db_releases['superfamily']
        def superfamily_dir = "${datadir}/${appsConfig.superfamily.dir}/${superfamily_release}"

        SUPERFAMILY(
            ch_seqs,
            "${superfamily_dir}/${appsConfig.superfamily.hmm}",
            "${superfamily_dir}/${appsConfig.superfamily.selfhits}",
            "${superfamily_dir}/${appsConfig.superfamily.cla}",
            "${superfamily_dir}/${appsConfig.superfamily.model}",
            "${superfamily_dir}/${appsConfig.superfamily.pdbj95d}"
        )

        results = results.mix(SUPERFAMILY.out)
    }

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}