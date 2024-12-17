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
    ch_seqs             // channel of tuples (index, fasta file, json file)
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results = Channel.empty()

    ch_seqs
        .map { index, fasta, json -> tuple( index, fasta ) }
        .set { ch_fasta }

    ch_seqs
        .map { index, fasta, json -> tuple( index, json ) }
        .set { ch_json }

    if (applications.contains("antifam")) {
        RUN_ANTIFAM(
            ch_fasta,
            "${datadir}/${appsConfig.antifam.hmm}"
        )

        PARSE_ANTIFAM(RUN_ANTIFAM.out)
        results = results.mix(PARSE_ANTIFAM.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        // Search Gene3D profiles
        SEARCH_GENE3D(
            ch_fasta,
            "${datadir}/${appsConfig.cathgene3d.hmm}")
        ch_gene3d = SEARCH_GENE3D.out

        // Select best domain matches
        RESOLVE_GENE3D(ch_gene3d)

        // Assign CATH superfamily to matches
        ASSIGN_CATH(
            RESOLVE_GENE3D.out,
            "${datadir}/${appsConfig.cathgene3d.model2sfs}",
            "${datadir}/${appsConfig.cathgene3d.disc_regs}"
        )

        // Join results and parse them
        ch_cathgene3d = ch_gene3d.join(ASSIGN_CATH.out)
        PARSE_CATHGENE3D(ch_cathgene3d)

        if (applications.contains("cathgene3d")) {
            results = results.mix(PARSE_CATHGENE3D.out)
        }

        if (applications.contains("cathfunfam")) {
            // Find unique CATH superfamilies with at least one hit
            PREPARE_FUNFAM(
                PARSE_CATHGENE3D.out,
                "${datadir}/${appsConfig.cathfunfam.dir}"
            )

            /*
                Join input fasta file with superfamilies.
                We split in smaller chunks to parallelize searching FunFam profiles
            */
            ch_fasta
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
                "${datadir}/${appsConfig.cathfunfam.dir}"
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
        RUN_RPSBLAST(
            ch_fasta,
            "${datadir}/${appsConfig.cdd.rpsblast_db}"
        )

        RUN_RPSPROC(
            RUN_RPSBLAST.out,
            "${datadir}/${appsConfig.cdd.rpsproc_db}"
        )

        PARSE_RPSPROC(RUN_RPSPROC.out)
        results = results.mix(PARSE_RPSPROC.out)
    }

    if (applications.contains("coils")) {
        RUN_COILS(ch_fasta)
        PARSE_COILS(RUN_COILS.out)
        results = results.mix(PARSE_COILS.out)
    }

    if (applications.contains("deeptmhmm")) {
        RUN_DEEPTMHMM(
            ch_fasta,
            appsConfig.deeptmhmm.dir
        )
        PARSE_DEEPTMHMM(RUN_DEEPTMHMM.out)
        results = results.mix(PARSE_DEEPTMHMM.out)
    }

    if (applications.contains("hamap")) {
        PREPROCESS_HAMAP(
            ch_fasta,
            "${datadir}/${appsConfig.hamap.hmm}"
        )

        PREPARE_HAMAP(
            PREPROCESS_HAMAP.out.join(ch_json),
            "${datadir}/${appsConfig.hamap.dir}"
        )

        RUN_HAMAP(
            PREPARE_HAMAP.out,
            "${datadir}/${appsConfig.hamap.dir}"
        )

        PARSE_HAMAP(RUN_HAMAP.out)
        results = results.mix(PARSE_HAMAP.out)
    }

    if (applications.contains("mobidblite")) {
        RUN_MOBIDBLITE(ch_fasta)

        PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out)
        results = results.mix(PARSE_MOBIDBLITE.out)
    }

    if (applications.contains("ncbifam")) {
        RUN_NCBIFAM(
            ch_fasta,
            "${datadir}/${appsConfig.ncbifam.hmm}"
        )

        PARSE_NCBIFAM(RUN_NCBIFAM.out)
        results = results.mix(PARSE_NCBIFAM.out)
    }

    if (applications.contains("panther")) {
        SEARCH_PANTHER(
            ch_fasta,
            "${datadir}/${appsConfig.panther.hmm}"
        )
        ch_panther = SEARCH_PANTHER.out

        PREPARE_TREEGRAFTER(
            ch_panther,
            "${datadir}/${appsConfig.panther.msf}"
        )

        RUN_TREEGRAFTER(
            PREPARE_TREEGRAFTER.out.fasta,
            "${datadir}/${appsConfig.panther.msf}"
        )

        PARSE_PANTHER(
            PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out)
        )
        results = results.mix(PARSE_PANTHER.out)
    }

    if (applications.contains("phobius")) {
        SEARCH_PHOBIUS(
            ch_fasta,
            appsConfig.phobius.dir
        )

        PARSE_PHOBIUS(SEARCH_PHOBIUS.out)
        results = results.mix(PARSE_PHOBIUS.out)
    }

    if (applications.contains("pfam")) {
        SEARCH_PFAM(
            ch_fasta,
            "${datadir}/${appsConfig.pfam.hmm}"
        )

        PARSE_PFAM(
            SEARCH_PFAM.out,
            "${datadir}/${appsConfig.pfam.seed}",
            "${datadir}/${appsConfig.pfam.clan}",
            "${datadir}/${appsConfig.pfam.dat}"
        )
        results = results.mix(PARSE_PFAM.out)
    }

    if (applications.contains("pirsf")) {
        RUN_PIRSF(ch_fasta,
            "${datadir}/${appsConfig.pirsf.hmm}"
        )
        PARSE_PIRSF(RUN_PIRSF.out,
            "${datadir}/${appsConfig.pirsf.dat}")

        results = results.mix(PARSE_PIRSF.out)
    }

    if (applications.contains("pirsr")) {
        RUN_PIRSR(ch_fasta,
            "${datadir}/${appsConfig.pirsr.hmm}")

        PARSE_PIRSR(RUN_PIRSR.out,
            "${datadir}/${appsConfig.pirsr.rules}")

        results = results.mix(PARSE_PIRSR.out)
    }

    if (applications.contains("prints")) {
        RUN_PRINTS(
            ch_fasta,
            "${datadir}/${appsConfig.prints.pval}"
        )

        PARSE_PRINTS(
            RUN_PRINTS.out,
            "${datadir}/${appsConfig.prints.hierarchy}"
        )
        results = results.mix(PARSE_PRINTS.out)
    }

    if (applications.contains("prositepatterns")) {
        RUN_PFSCAN(
            ch_fasta,
            "${datadir}/${appsConfig.prositepatterns.dat}",
            "${datadir}/${appsConfig.prositepatterns.evaluator}")

        PARSE_PFSCAN(RUN_PFSCAN.out)
        results = results.mix(PARSE_PFSCAN.out)
    }

    if (applications.contains("prositeprofiles")) {
        RUN_PFSEARCH(
            ch_fasta,
            "${datadir}/${appsConfig.prositeprofiles.dir}")
        PARSE_PFSEARCH(
            RUN_PFSEARCH.out,
            "${datadir}/${appsConfig.prositeprofiles.skip_flagged_profiles}")
        results = results.mix(PARSE_PFSEARCH.out)
    }

    if (applications.contains("sfld")) {
        RUN_SFLD(ch_fasta,
            "${datadir}/${appsConfig.sfld.hmm}")

        POST_PROCESS_SFLD(RUN_SFLD.out,
            "${datadir}/${appsConfig.sfld.sites_annotation}")

        PARSE_SFLD(POST_PROCESS_SFLD.out,
            "${datadir}/${appsConfig.sfld.hierarchy}")

        results = results.mix(PARSE_SFLD.out)
    }

    if (applications.contains("signalp_euk")) {
        RUN_SIGNALP_EUK(
            ch_fasta, 
            appsConfig.signalp_euk.organism, 
            appsConfig.signalp_euk.mode,
            appsConfig.signalp_euk.dir
        )
        PARSE_SIGNALP_EUK(RUN_SIGNALP_EUK.out)
        results = results.mix(PARSE_SIGNALP_EUK.out)
    }

    if (applications.contains("signalp_prok")) {
        RUN_SIGNALP_PROK(
            ch_fasta, 
            appsConfig.signalp_prok.organism, 
            appsConfig.signalp_prok.mode,
            appsConfig.signalp_prok.dir
        )
        PARSE_SIGNALP_PROK(RUN_SIGNALP_PROK.out)
        results = results.mix(PARSE_SIGNALP_PROK.out)
    }

    if (applications.contains("smart")) {
        SEARCH_SMART(
            ch_fasta,
            "${datadir}/${appsConfig.smart.hmmbin}"
        )

        PARSE_SMART(
            SEARCH_SMART.out.join(ch_json),
            "${datadir}/${appsConfig.smart.hmm}"
        )
        results = results.mix(PARSE_SMART.out)
    }

    if (applications.contains("superfamily")) {
        SEARCH_SUPERFAMILY(
            ch_fasta,
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
