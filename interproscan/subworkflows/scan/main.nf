include { RUN_ANTIFAM; PARSE_ANTIFAM                                              } from  "../../modules/antifam"
include { SEARCH_GENE3D; RESOLVE_GENE3D; ASSIGN_CATH_SUPFAM; PARSE_CATHGENE3D     } from  "../../modules/cath/gene3d"
include { PREPARE_FUNFAM     } from  "../../modules/cath/funfam"
include { RUN_RPSBLAST; RUN_RPSPROC; PARSE_RPSPROC                                } from  "../../modules/cdd"
include { RUN_COILS; PARSE_COILS                                                  } from  "../../modules/coils"
include { PREPROCESS_HAMAP; PREPARE_HAMAP; RUN_HAMAP; PARSE_HAMAP                 } from  "../../modules/hamap"
include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE                                        } from  "../../modules/mobidblite"
include { RUN_NCBIFAM; PARSE_NCBIFAM                                              } from  "../../modules/ncbifam"

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
        .map { index, fasta, json -> json }
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
        // SEARCH_GENE3D(
        //     ch_fasta,
        //     "${datadir}/${appsConfig.cathgene3d.hmm}")

        // SEARCH_GENE3D.out.view()

        // RESOLVE_GENE3D(SEARCH_GENE3D.out)

        // ASSIGN_CATH_SUPFAM(
        //     RESOLVE_GENE3D.out,
        //     "${datadir}/${appsConfig.cathgene3d.model2sfs}",
        //     "${datadir}/${appsConfig.cathgene3d.disc_regs}")

        // ch_cathgene3d = SEARCH_GENE3D.out.join(ASSIGN_CATH_SUPFAM.out)

        def dummy = channel.from(
            [tuple("1", file("/home/mblum/Projects/i6/work/a1/b38678716d6b10289423feeca79331/hmmsearch.out"))],
        )
        RESOLVE_GENE3D(dummy)
        ASSIGN_CATH_SUPFAM(
            RESOLVE_GENE3D.out,
            "${datadir}/${appsConfig.cathgene3d.model2sfs}",
            "${datadir}/${appsConfig.cathgene3d.disc_regs}")
        ch_cathgene3d = dummy.join(ASSIGN_CATH_SUPFAM.out)

        PARSE_CATHGENE3D(ch_cathgene3d)
        results = results.mix(PARSE_CATHGENE3D.out)

        if (applications.contains("cathfunfam")) {
            // TODO CATH-FunFam
            PREPARE_FUNFAM(PARSE_CATHGENE3D.out,
                "${datadir}/${appsConfig.cathfunfam.dir}")
        }
    }

    if (applications.contains("cdd")) {
        RUN_RPSBLAST(
            ch_fasta,
            "${datadir}/${appsConfig.cdd.rpsblast_db}")
        RUN_RPSPROC(
            RUN_RPSBLAST.out,
            "${datadir}/${appsConfig.cdd.rpsproc_db}")

        PARSE_RPSPROC(RUN_RPSPROC.out)
        results = results.mix(PARSE_RPSPROC.out)
    }

    if (applications.contains("coils")) {
        RUN_COILS(ch_fasta)
        PARSE_COILS(RUN_COILS.out)
        results = results.mix(PARSE_COILS.out)
    }

    if (applications.contains("hamap")) {
        PREPROCESS_HAMAP(ch_fasta,
            "${datadir}/${appsConfig.hamap.hmm}")

        PREPARE_HAMAP(PREPROCESS_HAMAP.out,
            ch_json,
            "${datadir}/${appsConfig.hamap.dir}")

        RUN_HAMAP(PREPARE_HAMAP.out,
            "${datadir}/${appsConfig.hamap.dir}")

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
        // TODO
    }

    if (applications.contains("phobius")) {
        // TODO
    }

    if (applications.contains("pfam")) {
        // TODO
    }

    if (applications.contains("pirsf")) {
        // TODO
    }

    if (applications.contains("pirsr")) {
        // TODO
    }

    if (applications.contains("prints")) {
        // TODO
    }

    if (applications.contains("prositepatterns")) {
        // TODO
    }

    if (applications.contains("prositeprofiles")) {
        // TODO
    }

    if (applications.contains("sfld")) {
        // TODO
    }

    if (applications.contains("smart")) {
        // TODO
    }

    if (applications.contains("superfamily")) {
        // TODO
    }

    if (applications.contains("signalp")) {
        // TODO
    }

    if (applications.contains("signalp_euk")) {
        // TODO
    }

    if (applications.contains("tmhmm")) {
        // TODO
    }

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}