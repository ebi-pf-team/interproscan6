include { RUN_PRINTS; PARSE_PRINTS } from  "../../../modules/prints"

workflow PRINTS {
    take:
    ch_seqs
    prints_pval
    prints_hierarchy

    main:
    RUN_PRINTS(
        ch_seqs,
        prints_pval
    )

    ch_prints = PARSE_PRINTS(
        RUN_PRINTS.out,
        prints_hierarchy
    )

    emit:
    ch_prints
}