include { SEARCH_PANTHER; PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "../../modules/panther"

workflow PANTHER {
    take:
    ch_fasta       // [meta, fasta]
    panther        // [hmm, msf]
    batch_size     // int, number of sequences per sub batch for searching

    main:
    ch_split = ch_fasta
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    ch_hmm = panther.map { hmm, _msf -> hmm }
    SEARCH_PANTHER(
        ch_split,
        ch_hmm
    )

    PREPARE_TREEGRAFTER(
        SEARCH_PANTHER.out,
    )

    ch_msf = panther.map { _hmm, msf -> msf }
    RUN_TREEGRAFTER(
        PREPARE_TREEGRAFTER.out,
        ch_msf
    )

    PARSE_PANTHER(
        PREPARE_TREEGRAFTER.out.join(
            RUN_TREEGRAFTER.out, 
            by: [0, 1],
            failOnDuplicate: true,
            failOnMismatch: true
        )
    )
        
    emit:
    PARSE_PANTHER.out
        .map { meta, _meta2, json -> tuple (meta, json) }  // [ meta, json ]
}