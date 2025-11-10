include { PREPARE_TMBED; RUN_TMBED_CPU; RUN_TMBED_GPU; PARSE_TMBED } from "../../modules/tmbed"

workflow TMBED {
    take:
    ch_seqs
    use_gpu
    subbatch_size       // params.subBatchSize
    chunk_size          // max number of aa per chunked sequence
    chunk_overlap       // overlap between two chunks of a sequence
    smooth_window       // length of the sliding window used during smoothing
    tmbed_batch_size    // max number of aa per TMbed batch

    main:
    PREPARE_TMBED(
        ch_seqs,
        chunk_size,
        chunk_overlap,
        use_gpu ? subbatch_size * 10 : subbatch_size.intdiv(5)
    )

    ch_split = PREPARE_TMBED.out
        .flatMap { meta, files ->
            files.collect { f -> tuple(meta, f) }
        }

    if (use_gpu) {
        RUN_TMBED_GPU(ch_split, tmbed_batch_size)
        ch_tmbed = RUN_TMBED_GPU.out
    } else {
        RUN_TMBED_CPU(ch_split, tmbed_batch_size)
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