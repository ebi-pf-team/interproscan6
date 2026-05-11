include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE } from  "../../modules/mobidblite"

workflow MOBIDBLITE {
    take:
    ch_seqs     // channel of tuples (meta, fasta file)
    batch_size  // int, number of sequences per sub batch for searching

    main:
    ch_split = ch_seqs
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    RUN_MOBIDBLITE(ch_split)

    PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out)

    emit:
    PARSE_MOBIDBLITE.out
        .map { meta, meta2, json -> tuple (meta, json) }  // [ meta, json ]
}
