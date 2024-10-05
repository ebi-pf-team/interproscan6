include { RUN_ANTIFAM; PARSE_ANTIFAM                                              } from  "../../modules/antifam"
include { SEARCH_GENE3D; RESOLVE_GENE3D; ASSIGN_CATH; PARSE_CATHGENE3D            } from  "../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; RESOLVE_FUNFAM; PARSE_FUNFAM             } from  "../../modules/cath/funfam"
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
            "${datadir}/${appsConfig.cathgene3d.disc_regs}")

        // Join results and parse them
        ch_cathgene3d = ch_gene3d.join(ASSIGN_CATH.out)
        PARSE_CATHGENE3D(ch_cathgene3d)

        if (applications.contains("cathgene3d")) {
            results = results.mix(PARSE_CATHGENE3D.out)
        }

        if (applications.contains("cathfunfam")) {
            // Find unique CATH superfamilies with at least one hit
            PREPARE_FUNFAM(PARSE_CATHGENE3D.out,
                "${datadir}/${appsConfig.cathfunfam.dir}")

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
            SEARCH_FUNFAM(ch_funfams,
                "${datadir}/${appsConfig.cathfunfam.dir}")

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