include { RUN_HMMER as SEARCH_PANTHER                         } from  "../../modules/hmmer"
include { PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "../../modules/panther"

workflow PANTHER {
    take:
    ch_seqs
    dir
    hmm
    msf
    
    main:
    SEARCH_PANTHER(
        ch_seqs,
        dir,
        hmm,
        "-Z 65000000 -E 0.001"
    )
    ch_panther = SEARCH_PANTHER.out

    PREPARE_TREEGRAFTER(
        ch_panther,
        dir,
        msf
    )

    RUN_TREEGRAFTER(
        PREPARE_TREEGRAFTER.out.fasta,
        dir,
        msf
    )

    ch_panther = PARSE_PANTHER(
        PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out)
    )
        
    emit:
    ch_panther
}