// codenarc-disable ModuleIncludedTwiceRule
include { RUN_SIGNALP_CPU as RUN_SIGNALP_CPU_EUK      } from  "../../modules/signalp"
include { RUN_SIGNALP_GPU as RUN_SIGNALP_GPU_EUK      } from  "../../modules/signalp"
include { PARSE_SIGNALP as PARSE_SIGNALP_EUK          } from  "../../modules/signalp"
include { RUN_SIGNALP_CPU as RUN_SIGNALP_CPU_PROK     } from  "../../modules/signalp"
include { RUN_SIGNALP_GPU as RUN_SIGNALP_GPU_PROK     } from  "../../modules/signalp"
include { PARSE_SIGNALP as PARSE_SIGNALP_PROK         } from  "../../modules/signalp"

workflow SIGNALP {
    take:
    ch_seqs        // channel of tuples (index, fasta file)
    applications   // list of applications to run: signalp_euk, signalp_prok
    euk_organism   // str, signalp organism option, always 'eukarya' for eukaryotes
    euk_mode       // str, signalp mode option, e.g. 'fast' or 'slow'
    euk_dir        // str repr of the path to the eukaryotic signalp data dir
    euk_gpu        // boolean, use GPU for eukaryotic signalp
    prok_organism  // str, signalp organism option, always 'prokarya' for prokaryotes
    prok_mode      // str, signalp mode option, e.g. 'fast' or 'slow'
    prok_dir       // str repr of the path to the prokaryotic signalp data dir
    prok_gpu       // boolean, use GPU for prokaryotic signalp
    batch_size     // int, number of sequences per sub-batch for searching

    main:
    results = Channel.empty()

    if (applications.contains("signalp_euk")) {
        if (euk_gpu) {
            RUN_SIGNALP_GPU_EUK(
                ch_seqs,
                euk_organism,
                euk_mode,
                euk_dir
            )
            ch_euk = RUN_SIGNALP_GPU_EUK.out
        } else {
            ch_euk_split = ch_seqs
                .splitFasta( by: batch_size * 10, file: true )
            RUN_SIGNALP_CPU_EUK(
                ch_euk_split,
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
        if (prok_gpu) {
            RUN_SIGNALP_GPU_PROK(
                ch_seqs,
                prok_organism,
                prok_mode,
                prok_dir
            )
            ch_prok = RUN_SIGNALP_GPU_PROK.out
        } else {
            ch_prok_split = ch_seqs
                .splitFasta( by: batch_size * 10, file: true )
            RUN_SIGNALP_CPU_PROK(
                ch_prok_split,
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