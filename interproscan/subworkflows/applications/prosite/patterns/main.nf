include { RUN_PFSCAN; PARSE_PFSCAN } from  "../../../../modules/prosite/patterns"

workflow PROSITE_PATTERNS {
    take:
    ch_seqs
    prositepatterns_dat
    prositepatterns_evaluator

    main:
    RUN_PFSCAN(
        ch_seqs,
        prositepatterns_dat,
        prositepatterns_evaluator
    )

    ch_prosite = PARSE_PFSCAN(RUN_PFSCAN.out)

    emit:
    ch_prosite
}