include { RUN_DEEPTMHMM_CPU; RUN_DEEPTMHMM_GPU; PARSE_DEEPTMHMM } from  "../../modules/deeptmhmm"

workflow DEEPTMHMM {
    take:
    ch_seqs       // channel of tuples (index, fasta file)
    deeptmhmm_dir // str repr of the data directory path
    use_gpu       // boolean to run on GPU
    batch_size    // int, number of sequences per sub batch for searching

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
