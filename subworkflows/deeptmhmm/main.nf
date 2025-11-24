include { RUN_DEEPTMHMM_CPU; RUN_DEEPTMHMM_GPU; PARSE_DEEPTMHMM } from  "../../modules/deeptmhmm"

workflow DEEPTMHMM {
    take:
    ch_seqs
    deeptmhmm_dir
    use_gpu
    batch_size

    main:
    if (use_gpu) {
        ch_split = ch_seqs
            .splitFasta( by: batch_size * 2, file: true )
        RUN_DEEPTMHMM_GPU(
            ch_split,
            deeptmhmm_dir
        )
        ch_deeptmhmm = RUN_DEEPTMHMM_GPU.out
    } else {
        ch_split = ch_seqs
            .splitFasta( by: batch_size.intdiv(5), file: true )

        RUN_DEEPTMHMM_CPU(
            ch_split,
            deeptmhmm_dir
        )
        ch_deeptmhmm = RUN_DEEPTMHMM_CPU.out
    }
    
    ch_results = PARSE_DEEPTMHMM(ch_deeptmhmm)

    emit:
    ch_results
}
