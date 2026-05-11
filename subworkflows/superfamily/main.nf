include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY } from  "../../modules/superfamily"

workflow SUPERFAMILY {
    take:
    fasta      // [meta, fasta]
    ssf        // [dir, hmm, selfhits, cla, model, pdbj95d]
    batch_size // [int]

    main:
    ch_split = fasta
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    SEARCH_SUPERFAMILY(ch_split, ssf)

    ch_parse = ssf.map { dir, hmm, selfhits, cla, model, pdbj95d -> tuple (dir, hmm, model) }

    PARSE_SUPERFAMILY(
        SEARCH_SUPERFAMILY.out,
        ch_parse
    )

    emit:
    PARSE_SUPERFAMILY.out
        .map { meta, meta2, json -> tuple (meta, json) }  // [ meta, json ]
}