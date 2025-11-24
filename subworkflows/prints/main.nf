include { RUN_PRINTS; PARSE_PRINTS } from  "../../modules/prints"

workflow PRINTS {
    take:
    ch_seqs             // channel of tuples (index, fasta file)
    dirpath             // data directory path
    pvalfile            // path to the p-value file in the data dir -> datadir/pvalFile
    hierarchyfile       // path to the hierarchy file in the data dir -> datadir/hierarchyFile
    batch_size          // number of sequences per sub-batch for searching

    main:
    ch_split = ch_seqs
        .splitFasta( by: batch_size, file: true )

    RUN_PRINTS(
        ch_split,
        dirpath,
        pvalfile
    )

    ch_prints = PARSE_PRINTS(
        RUN_PRINTS.out,
        dirpath,
        hierarchyfile
    )

    emit:
    ch_prints
}