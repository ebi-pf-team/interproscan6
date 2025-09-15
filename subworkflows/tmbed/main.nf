include { RUN_TMBED_CPU; RUN_TMBED_GPU; PARSE_TMBED } from "../../modules/tmbed"

workflow TMBED {
    take:
    ch_seqs
    use_gpu
    batch_size

    main:
    if (use_gpu) {
        RUN_TMBED_GPU(
            ch_seqs
        )
        ch_tmbed = RUN_TMBED_GPU.out
    } else {
        ch_split = ch_seqs
            .splitFasta( by: batch_size, file: true )

        RUN_TMBED_CPU (
            ch_split
        )
        ch_tmbed = RUN_TMBED_CPU.out
    }

    ch_results = PARSE_TMBED(ch_tmbed)

    emit:
    ch_results
}