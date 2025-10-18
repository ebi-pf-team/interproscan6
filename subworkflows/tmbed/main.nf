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
    PREPARE_TMBED(
        ch_seqs,
        chunk_size,
        chunk_overlap,
        use_gpu ? 5000 : 100 // max number of sequences per FASTA file
    )

    ch_split = PREPARE_TMBED.out
        .flatMap { meta, files ->
            files.collect { f -> tuple(meta, f) }
        }

    if (use_gpu) {
        RUN_TMBED_GPU(ch_split, batch_size)
        ch_tmbed = RUN_TMBED_GPU.out
    } else {
        RUN_TMBED_CPU(ch_split, batch_size)
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