include { RUN_COILS; PARSE_COILS } from  "${moduleDir}/../../modules/coils"

workflow COILS {
    take:
    ch_seqs // channel of tuples (index, fasta file)

    main:
    RUN_COILS(ch_seqs)
    ch_coils = PARSE_COILS(RUN_COILS.out)

    emit:
    ch_coils
}