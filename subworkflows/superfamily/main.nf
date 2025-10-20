include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY } from  "../../modules/superfamily"

workflow SUPERFAMILY {
    take:
    ch_seqs
    dirpath
    hmm
    selfhits
    cla
    model
    pdbj95d
    batch_size

    main:
    results = Channel.empty()
    ch_split = ch_seqs
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    SEARCH_SUPERFAMILY(
        ch_split,
        dirpath,
        hmm,
        selfhits,
        cla,
        model,
        pdbj95d
    )

    results = results.mix(PARSE_SUPERFAMILY(
        SEARCH_SUPERFAMILY.out,
        dirpath,
        model,
        hmm
    ))

    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}