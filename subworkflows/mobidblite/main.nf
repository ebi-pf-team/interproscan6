include { RUN_MOBIDBLITE; PARSE_MOBIDBLITE } from  "../../modules/mobidblite"

workflow MOBIDBLITE {
    take:
    ch_seqs

    main:
    RUN_MOBIDBLITE(ch_seqs)

    ch_mobidblite =PARSE_MOBIDBLITE(RUN_MOBIDBLITE.out)

    emit:
    ch_mobidblite
}