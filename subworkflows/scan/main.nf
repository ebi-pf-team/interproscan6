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
    fasta            // [meta, fasta]
    appl_config
    appl_dirs  
    applications     // application to run in this workflow
    all_applications // all applications to run across all instances of SCAN_SEQUENCES
    batch_size       // int, sub-batch size for computationally demanding apps

    main:
    results = channel.empty()

    if (applications.contains("antifam")) {
        ch_antifam = appl_dirs
            .filter { name, dirpath -> name == "antifam" }
            .first()
            .map    { name, dirpath -> dirpath.resolve(appl_config.antifam.hmm) }

        ANTIFAM(fasta, ch_antifam)
        results = results.mix(ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        ch_cathgene3d = appl_dirs
            .filter { name, dirpath -> name == "cathgene3d" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.cathgene3d.hmm),
                    dirpath.resolve(appl_config.cathgene3d.model2sfs),
                    dirpath.resolve(appl_config.cathgene3d.disc_regs)
                )
            }

        ch_cathfunfam = appl_dirs
            .filter { name, dirpath -> name == "cathfunfam" }
            .first()
            .map    { name, dirpath -> dirpath }
        
        CATH(
            fasta,
            applications.contains("cathgene3d"),
            ch_cathgene3d,
            applications.contains("cathfunfam"),
            ch_cathfunfam,
            batch_size
        )

        results = results.mix(CATH.out)
    }

    if (applications.contains("cdd")) {
        ch_cdd = appl_dirs
            .filter { name, dirpath -> name == "cdd" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath,
                    appl_config.cdd.rpsblast_db,
                    appl_config.cdd.rpsproc_db
                )
            }

        CDD(fasta, ch_cdd)
        results = results.mix(CDD.out)
    }

    if (applications.contains("coils")) {
        COILS(fasta)
        results = results.mix(COILS.out)
    }

    if (applications.contains("deeptmhmm")) {
        DEEPTMHMM(
            fasta,
            appl_config.deeptmhmm.dir,
            appl_config.deeptmhmm.use_gpu,
            batch_size
        )
        results = results.mix(DEEPTMHMM.out)
    }

    if (applications.contains("hamap")) {
        ch_hamap = appl_dirs
            .filter { name, dirpath -> name == "hamap" }
            .first()
            .map    { name, dirpath -> dirpath.resolve(appl_config.hamap.profiles) }

        HAMAP(fasta, ch_hamap)
        results = results.mix(HAMAP.out)
    }

    if (applications.contains("interpro_n")) {
        INTERPRO_N(
            fasta,
            all_applications,
            appl_config.interpro_n.dir,
            appl_config.interpro_n.use_gpu,
            appl_config.interpro_n.batch_size,
            batch_size
        )

        results = results.mix(INTERPRO_N.out)
    }

    if (applications.contains("mobidblite")) {
        MOBIDBLITE(fasta, batch_size)
        results = results.mix(MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        ch_ncbifam = appl_dirs
            .filter { name, dirpath -> name == "ncbifam" }
            .first()
            .map    { name, dirpath -> dirpath.resolve(appl_config.ncbifam.hmm) }

        NCBIFAM(fasta, ch_ncbifam)
        results = results.mix(NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        ch_panther = appl_dirs
            .filter { name, dirpath -> name == "panther" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.panther.hmm),
                    dirpath.resolve(appl_config.panther.msf)
                )
            }

        PANTHER(fasta, ch_panther, batch_size)
        results = results.mix(PANTHER.out)
    }

    if (applications.contains("pfam")) {
        ch_pfam = appl_dirs
            .filter { name, dirpath -> name == "pfam" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.pfam.hmm),
                    dirpath.resolve(appl_config.pfam.dat)
                )
            }
        
        PFAM(fasta, ch_pfam)
        results = results.mix(PFAM.out)
    }

    if (applications.contains("phobius")) {
        PHOBIUS(
            fasta,
            appl_config.phobius.dir
        )
        results = results.mix(PHOBIUS.out)
    }

    if (applications.contains("pirsf")) {
        ch_pirsf = appl_dirs
            .filter { name, dirpath -> name == "pirsf" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.pirsf.hmm),
                    dirpath.resolve(appl_config.pirsf.dat)
                )
            }

        PIRSF(fasta, ch_pirsf)
        results = results.mix(PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        ch_pirsr = appl_dirs
            .filter { name, dirpath -> name == "pirsr" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.pirsr.hmm),
                    dirpath.resolve(appl_config.pirsr.rules)
                )
            }

        PIRSR(fasta, ch_pirsr)
        results = results.mix(PIRSR.out)
    }

    if (applications.contains("prints")) {
        ch_prints = appl_dirs
            .filter { name, dirpath -> name == "prints" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.prints.pval),
                    dirpath.resolve(appl_config.prints.hierarchy)
                )
            }

        PRINTS(fasta, ch_prints, batch_size)
        results = results.mix(PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        ch_prositepatterns = appl_dirs
            .filter { name, dirpath -> name == "prositepatterns" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.prositepatterns.dat),
                    dirpath.resolve(appl_config.prositepatterns.evaluator)
                )
            }

        PROSITE_PATTERNS(fasta, ch_prositepatterns)
        results = results.mix(PROSITE_PATTERNS.out)
    }

    if (applications.contains("prositeprofiles")) {
        ch_prositeprofiles = appl_dirs
            .filter { name, dirpath -> name == "prositeprofiles" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.prositeprofiles.profiles),
                    dirpath.resolve(appl_config.prositeprofiles.skip_flagged_profiles)
                )
            }

        PROSITE_PROFILES(fasta, ch_prositeprofiles)
        results = results.mix(PROSITE_PROFILES.out)
    }

    if (applications.contains("sfld")) {
        ch_sfld = appl_dirs
            .filter { name, dirpath -> name == "sfld" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.sfld.hmm),
                    dirpath.resolve(appl_config.sfld.sites_annotation),
                    dirpath.resolve(appl_config.sfld.hierarchy)
                )
            }
    
        SFLD(fasta, ch_sfld)
        results = results.mix(SFLD.out)
    }

    if (applications.contains("signalp_euk") || applications.contains("signalp_prok")) {
        SIGNALP(
            fasta,
            applications,
            appl_config.signalp_euk.organism,
            appl_config.signalp_euk.mode,
            appl_config.signalp_euk.dir,
            appl_config.signalp_euk.use_gpu,
            appl_config.signalp_prok.organism,
            appl_config.signalp_prok.mode,
            appl_config.signalp_prok.dir,
            appl_config.signalp_prok.use_gpu,
            batch_size
        )
        results = results.mix(SIGNALP.out)
    }

    if (applications.contains("smart")) {
        ch_smart = appl_dirs
            .filter { name, dirpath -> name == "smart" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath.resolve(appl_config.smart.hmm3),
                    dirpath.resolve(appl_config.smart.hmm2)
                )
            }

        SMART(fasta, ch_smart)
        results = results.mix(SMART.out)
    }

    if (applications.contains("superfamily")) {
        ch_superfamily = appl_dirs
            .filter { name, dirpath -> name == "superfamily" }
            .first()
            .map    { name, dirpath ->
                tuple(
                    dirpath,
                    appl_config.superfamily.hmm,
                    appl_config.superfamily.selfhits,
                    appl_config.superfamily.cla,
                    appl_config.superfamily.model,
                    appl_config.superfamily.pdbj95d,
                )
            }

        SUPERFAMILY(fasta, ch_superfamily, batch_size)
        results = results.mix(SUPERFAMILY.out)
    }

    if (applications.contains("tmbed")) {
        TMBED(
            fasta,
            appl_config.tmbed.use_gpu,
            appl_config.tmbed.chunk_size,
            appl_config.tmbed.chunk_overlap,
            appl_config.tmbed.smooth_window,
            appl_config.tmbed.batch_size,
            batch_size
        )

        results = results.mix(TMBED.out)
    }

    // Add a dummy null value for each fasta so there are no mismatches when calling join()
    results = results.mix(fasta.map { meta, fasta -> tuple(meta, null) })

    results = fasta
        .join(
            results.groupTuple(),
            failOnDuplicate: true,
            failOnMismatch: true
        )
        .map { meta, fasta, jsons ->
            tuple(meta, fasta, jsons.collect().findAll { f -> f != null })
        }

    no_matches = REPORT_NO_MATCHES(results)

    results = results.join(
        no_matches,
        failOnDuplicate: true,
        failOnMismatch: true
    )
    .map { meta, fasta, jsons, json ->
        tuple(meta, jsons + [json])
    }

    emit:
    results
}