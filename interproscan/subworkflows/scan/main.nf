include { RUN_ANTIFAM           } from  "../../modules/antifam"
include { PARSE_ANTIFAM         } from  "../../modules/antifam"
include { RUN_MOBIDBLITE        } from  "../../modules/mobidblite"
include { PARSE_MOBIDBLITE      } from  "../../modules/mobidblite"

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

    results
        .groupTuple()
        .set { grouped_results }

    emit:
    grouped_results
}