include { RUN_PRINTS; PARSE_PRINTS } from  "../../modules/prints"

workflow PRINTS {
    take:
    fasta       // [meta, fasta]
    prints      // [pval, hierarchy]
    batch_size  // [int]

    main:
    ch_split = fasta
        .splitFasta( by: batch_size, file: true )

    pval = prints.map { pval, hierarchy -> pval }
    hierarchy = prints.map { pval, hierarchy -> hierarchy }

    RUN_PRINTS(
        ch_split,
        pval
    )

    PARSE_PRINTS(
        RUN_PRINTS.out,
        hierarchy
    )

    emit:
    PARSE_PRINTS.out
}