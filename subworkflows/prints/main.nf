include { RUN_PRINTS; PARSE_PRINTS } from  "../../modules/prints"

workflow PRINTS {
    take:
    ch_seqs
    dirpath
    pvalfile
    hierarchyfile
    batch_size

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