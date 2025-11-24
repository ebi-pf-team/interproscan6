include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE } from  "../../modules/mobidblite"

workflow MOBIDBLITE {
    take:
    ch_seqs     // channel of tuples (meta, fasta file)
    batch_size  // number of sequences per sub batch for searching

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

    RUN_MOBIDBLITE(ch_split)

    results = results.mix(PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out))

    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}
