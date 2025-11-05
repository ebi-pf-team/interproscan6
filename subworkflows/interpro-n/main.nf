include { PREPARE_INTERPRO_N; RUN_INTERPRO_N_CPU; RUN_INTERPRO_N_GPU; PARSE_INTERPRO_N } from  "../../modules/interpro-n"

workflow INTERPRO_N {
    take:
    ch_seqs
    applications
    checkpoint_dir
    use_gpu
    res_batch_size      // Max number of residues per InterPro-N batch
    seq_batch_size      // params.subBatchSize

    main:
    PREPARE_INTERPRO_N(
        ch_seqs, 
        use_gpu ? seq_batch_size * 10 : seq_batch_size.intdiv(5)
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

    ch_results = PARSE_INTERPRO_N(ch_json, applications)

    emit:
    ch_results
}