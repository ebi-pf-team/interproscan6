include { RUN_DEEPTMHMM_CPU; RUN_DEEPTMHMM_GPU; PARSE_DEEPTMHMM } from  "../../modules/deeptmhmm"

workflow DEEPTMHMM {
    take:
    ch_seqs
    deeptmhmm_dir
    use_gpu

    main:
    if (use_gpu) {
        RUN_DEEPTMHMM_GPU(
            ch_seqs,
            deeptmhmm_dir
        )
        ch_deeptmhmm = RUN_DEEPTMHMM_GPU.out
    } else {
        ch_split = ch_seqs
            .splitFasta( by: 100, file: true )

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