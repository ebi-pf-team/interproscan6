include { PREPARE_INTERPRO_N; RUN_INTERPRO_N_CPU; RUN_INTERPRO_N_GPU; PARSE_INTERPRO_N } from  "../../modules/interpro_n"

workflow INTERPRO_N {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    applications        // list of applications (str) to report InterPro-N matches
    checkpoint_dir      // str repr of the path to the checkpoint directory
    use_gpu             // str repr of the boolean to run on GPU
    res_batch_size      // int, max number of residues per InterPro-N batch
    seq_batch_size      // int, number of sequences per sub batch for searching

    main:
    def max_length    = 2047  // Maximum sequence length
    def chunk_overlap = 1000  // Number of overlapping residues between consecutive chunks
    def match_overlap = 0.25  // Fractional overlap threshold for merging matches
    PREPARE_INTERPRO_N(
        ch_seqs, 
        use_gpu ? seq_batch_size * 10 : seq_batch_size.intdiv(5),
        max_length,
        chunk_overlap
    )

    ch_tsv = PREPARE_INTERPRO_N.out
        .flatMap { meta, files ->
            files.collect { f -> tuple(meta, f) }
        }

    if (use_gpu) {
        RUN_INTERPRO_N_GPU(
            ch_tsv,
            checkpoint_dir,
            res_batch_size
        )

        ch_json = RUN_INTERPRO_N_GPU.out
    } else {
        RUN_INTERPRO_N_CPU(
            ch_tsv,
            checkpoint_dir,
            res_batch_size
        )
        ch_json = RUN_INTERPRO_N_CPU.out
    }

    // Pair each batch of results with its input sequences, whose lengths are
    // needed to bound match coordinates
    PARSE_INTERPRO_N(
        ch_json.combine(ch_seqs, by: 0),
        applications,
        max_length,
        chunk_overlap,
        match_overlap)

    emit:
    PARSE_INTERPRO_N.out  // [ meta, json ]
}