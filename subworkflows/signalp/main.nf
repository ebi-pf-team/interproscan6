include { RUN_SIGNALP_CPU as RUN_SIGNALP_CPU_EUK      } from  "../../modules/signalp"
include { RUN_SIGNALP_GPU as RUN_SIGNALP_GPU_EUK      } from  "../../modules/signalp"
include { PARSE_SIGNALP as PARSE_SIGNALP_EUK          } from  "../../modules/signalp"
include { RUN_SIGNALP_CPU as RUN_SIGNALP_CPU_PROK     } from  "../../modules/signalp"
include { RUN_SIGNALP_GPU as RUN_SIGNALP_GPU_PROK     } from  "../../modules/signalp"
include { PARSE_SIGNALP as PARSE_SIGNALP_PROK         } from  "../../modules/signalp"

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
    use_gpu

    main:
    results = Channel.empty()

    if (applications.contains("signalp_euk")) {
        if (use) {
            RUN_SIGNALP_GPU_EUK(
                ch_seqs,
                euk_organism,
                euk_mode,
                euk_dir
            )
            ch_euk = RUN_SIGNALP_GPU_EUK.out
        } else {
            RUN_SIGNALP_CPU_EUK(
                ch_seqs,
                euk_organism,
                euk_mode,
                euk_dir
            )
            ch_euk = RUN_SIGNALP_CPU_EUK.out
        }
        PARSE_SIGNALP_EUK(ch_euk)
        results = results.mix(PARSE_SIGNALP_EUK.out)
    }

    if (applications.contains("signalp_prok")) {
        if (use_gpu) {
            RUN_SIGNALP_GPU_PROK(
                ch_seqs,
                prok_organism,
                prok_mode,
                prok_dir
            )
            ch_prok = RUN_SIGNALP_GPU_PROK.out
        } else {
            RUN_SIGNALP_CPU_PROK(
                ch_seqs,
                prok_organism,
                prok_mode,
                prok_dir
            )
            ch_prok = RUN_SIGNALP_CPU_PROK.out
        }
        PARSE_SIGNALP_PROK(ch_prok)
        results = results.mix(PARSE_SIGNALP_PROK.out)
    }

    emit:
    results
}