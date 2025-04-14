include { RUN_SIGNALP as RUN_SIGNALP_EUK; PARSE_SIGNALP as PARSE_SIGNALP_EUK      } from  "../../../modules/signalp"
include { RUN_SIGNALP as RUN_SIGNALP_PROK; PARSE_SIGNALP as PARSE_SIGNALP_PROK    } from  "../../../modules/signalp"

workflow SIGNALP {
    take:
    ch_seqs
    applications
    euk_organism
    euk_mode
    euk_dir
    prok_organism
    prok_mode
    prok_dir

    main:
    results = Channel.empty()

    if (applications.contains("signalp_euk")) {
        RUN_SIGNALP_EUK(
            ch_seqs,
            euk_organism,
            euk_mode,
            euk_dir
        )
        PARSE_SIGNALP_EUK(RUN_SIGNALP_EUK.out)
        results = results.mix(PARSE_SIGNALP_EUK.out)
    }

    if (applications.contains("signalp_prok")) {
        RUN_SIGNALP_PROK(
            ch_seqs,
            prok_organism,
            prok_mode,
            prok_dir
        )
        PARSE_SIGNALP_PROK(RUN_SIGNALP_PROK.out)
        results = results.mix(PARSE_SIGNALP_PROK.out)
    }

    emit:
    results
}