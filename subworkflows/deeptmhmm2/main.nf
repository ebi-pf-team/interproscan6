include { RUN_DEEPTMHMM2_CPU; RUN_DEEPTMHMM2_GPU; PARSE_DEEPTMHMM2 } from "../../modules/deeptmhmm2"

workflow DEEPTMHMM2 {
    take:
    ch_seqs       // channel of tuples (index, fasta file)
    use_gpu       // boolean to run on GPU
    batch_size    // int, number of sequences per sub batch for searching

    main:
    log.warn "DeepTMHMM2 is free for academic use. Commercial users should contact https://dtu.biolib.com/DeepTMHMM2 for licensing."
    
    if (use_gpu) {
        ch_split = ch_seqs
            .splitFasta( by: batch_size * 2, file: true )

        RUN_DEEPTMHMM2_GPU(ch_split)
        ch_deeptmhmm2 = RUN_DEEPTMHMM2_GPU.out
    } else {
        ch_split = ch_seqs
            .splitFasta( by: batch_size.intdiv(5), file: true )

        RUN_DEEPTMHMM2_CPU(ch_split)
        ch_deeptmhmm2 = RUN_DEEPTMHMM2_CPU.out
    }

    PARSE_DEEPTMHMM2(ch_deeptmhmm2)

    emit:
    PARSE_DEEPTMHMM2.out  // [ meta, json ]
}
