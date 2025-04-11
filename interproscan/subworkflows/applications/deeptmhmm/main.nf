include { RUN_DEEPTMHMM; PARSE_DEEPTMHMM } from  "../../../modules/deeptmhmm"

workflow DEEPTMHMM {
    take:
    ch_seqs
    deeptmhmm_dir

    main:
    RUN_DEEPTMHMM(
        ch_seqs,
        deeptmhmm_dir
    )
    ch_deeptmhmm = PARSE_DEEPTMHMM(RUN_DEEPTMHMM.out)

    emit:
    ch_deeptmhmm
}