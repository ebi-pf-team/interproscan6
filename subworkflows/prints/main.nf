include { RUN_PRINTS; PARSE_PRINTS } from  "../../modules/prints"

workflow PRINTS {
    take:
    ch_seqs
    dirpath
    pvalfile
    hierarchyfile

    main:
    ch_split = ch_seqs
        .splitFasta( by: 1000, file: true )

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