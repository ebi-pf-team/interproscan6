include { SEARCH_SUPERFAMILY; PARSE_SUPERFAMILY } from  "${moduleDir}/../../modules/superfamily"

workflow SUPERFAMILY {
    take:
    ch_seqs      // channel of tuples (meta, fasta file)
    dirpath      // str repr of the data directory path
    hmm          // str repr of the path to the HMM file in the data dir -> datadir/hmm
    selfhits     // str repr of the path to the selfhits file in the data dir -> datadir/selfhits
    cla          // str repr of the path to the cla file in the data dir -> datadir/cla
    model        // str repr of the path to the model dir in the data dir -> datadir/model
    pdbj95d      // str repr of the path to the pdbj95d file in the data dir -> datadir/pdbj95d
    batch_size   // int, number of sequences per sub-batch for searching

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