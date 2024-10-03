include { RUN_ANTIFAM; PARSE_ANTIFAM              } from  "../../modules/antifam"
include { RUN_COILS; PARSE_COILS                  } from  "../../modules/coils"
include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE        } from  "../../modules/mobidblite"

workflow SCAN_SEQUENCES {
    take:
    ch_fasta            // channel of tuples (index, fasta)
    applications        // list of applications to run
    appsConfig          // map of applications
    datadir             // path to data directory

    main:
    results          = Channel.empty()

    if (applications.contains("antifam")) {
        ch_out = RUN_ANTIFAM(
            ch_fasta,
            "${datadir}/${appsConfig.antifam.hmm}"
        )

        PARSE_ANTIFAM(ch_out)
        results = results.mix(PARSE_ANTIFAM.out)
    }

    if (applications.contains("mobidblite")) {
        ch_out = RUN_MOBIDBLITE(ch_fasta)
        PARSE_MOBIDBLITE(ch_out)
        results = results.mix(PARSE_MOBIDBLITE.out)
    }

    if (applications.contains("coils")) {
        ch_out = RUN_COILS(ch_fasta)
        PARSE_COILS(ch_out)
        results = results.mix(PARSE_COILS.out)
    }

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}