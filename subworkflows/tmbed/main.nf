include { PREPARE_TMBED; RUN_TMBED_CPU; RUN_TMBED_GPU; PARSE_TMBED } from "../../modules/tmbed"

workflow TMBED {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    use_gpu             // boolean to use GPU for prediction
    chunk_size          // int, max number of aa per chunked sequence
    chunk_overlap       // int, overlap between two chunks of a sequence
    smooth_window       // int, length of the sliding window used during smoothing
    tmbed_batch_size    // int, max number of aa per TMbed batch
    batch_size          // int, number of sequences per sub-batch for preparing

    main:
    PREPARE_TMBED(
        ch_seqs,
        chunk_size,
        chunk_overlap,
        use_gpu ? batch_size * 10 : batch_size.intdiv(5)
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

    PARSE_TMBED(
        ch_tmbed, 
        chunk_overlap, 
        smooth_window
    )

    emit:
    PARSE_TMBED.out  // [ meta, json ]
}