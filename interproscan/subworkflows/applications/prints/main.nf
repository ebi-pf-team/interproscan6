include { RUN_PRINTS; PARSE_PRINTS } from  "../../../modules/prints"

workflow PRINTS {
    take:
    ch_seqs
    dirpath
    pvalfile
    hierarchyfile

    main:
    RUN_PRINTS(
        ch_seqs,
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