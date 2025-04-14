include { PREPROCESS_HAMAP; PREPARE_HAMAP; RUN_HAMAP; PARSE_HAMAP } from  "../../../modules/hamap"

workflow HAMAP {
    take:
    ch_seqs
    hamap_hmm
    hamap_dir

    main:
    PREPROCESS_HAMAP(
        ch_seqs,
        hamap_hmm
    )

    PREPARE_HAMAP(
        PREPROCESS_HAMAP.out,
        hamap_dir
    )

    RUN_HAMAP(
        PREPARE_HAMAP.out,
        hamap_dir
    )

    ch_hamap = PARSE_HAMAP(RUN_HAMAP.out)

    emit:
    ch_hamap
}