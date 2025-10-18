include { PREPARE_TMBED; RUN_TMBED_CPU; RUN_TMBED_GPU; PARSE_TMBED } from "../../modules/tmbed"

workflow TMBED {
    take:
    ch_seqs
    use_gpu
    chunk_size
    chunk_overlap
    smooth_window
    batch_size

    main:
    if (use_gpu) {
        PREPARE_TMBED(
            ch_seqs,
            chunk_size,
            chunk_overlap
        )

        RUN_TMBED_GPU(PREPARE_TMBED.out, batch_size)
        ch_tmbed = RUN_TMBED_GPU.out
    } else {
        ch_split = ch_seqs
            .splitFasta( by: 100, file: true )

        PREPARE_TMBED(
            ch_split,
            chunk_size,
            chunk_overlap
        )

        RUN_TMBED_CPU(PREPARE_TMBED.out, batch_size)
        ch_tmbed = RUN_TMBED_CPU.out
    }

    ch_results = PARSE_TMBED(
        ch_tmbed, 
        chunk_overlap, 
        smooth_window
    )

    emit:
    ch_results
}