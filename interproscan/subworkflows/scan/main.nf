include { RUN_ANTIFAM; PARSE_ANTIFAM                  } from  "../../modules/antifam"
include { RUN_RPSBLAST; RUN_RPSPROC; PARSE_RPSPROC    } from  "../../modules/cdd"
include { RUN_COILS; PARSE_COILS                      } from  "../../modules/coils"
include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE            } from  "../../modules/mobidblite"

workflow SCAN_SEQUENCES {
    take:
    ch_fasta            // channel of tuples (index, fasta)
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results = Channel.empty()

    if (applications.contains("antifam")) {
        RUN_ANTIFAM(
            ch_fasta,
            "${datadir}/${appsConfig.antifam.hmm}"
        )

        PARSE_ANTIFAM(RUN_ANTIFAM.out)
        results = results.mix(PARSE_ANTIFAM.out)
    }

    if (applications.contains("mobidblite")) {
        RUN_MOBIDBLITE(ch_fasta)
        PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out)
        results = results.mix(PARSE_MOBIDBLITE.out)
    }

    if (applications.contains("cathgene3d") || applications.contains("cathfunfam")) {
        // TODO CATH-GENE3D

        if (applications.contains("cathfunfam")) {
            // TODO CATH-FunFam
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
        // TODO
    }

    if (applications.contains("ncbifam")) {
        // TODO
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