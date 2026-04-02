include { RUN_COILS; PARSE_COILS } from  "../../modules/coils"

workflow COILS {
    take:
    ch_seqs // channel of tuples (index, fasta file)

    main:
    RUN_COILS(ch_seqs)
    PARSE_COILS(RUN_COILS.out)

    emit:
    PARSE_COILS.out // [ meta, json ]
}