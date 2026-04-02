include { SEARCH_PANTHER; PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "../../modules/panther"

workflow PANTHER {
    take:
    fasta          // [meta, fasta]
    panther        // [hmm, msf]
    batch_size     // int, number of sequences per sub batch for searching

    main:
    ch_split = fasta
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    hmm = panther.map { hmm, msf -> hmm }
    SEARCH_PANTHER(
        ch_split,
        hmm
    )

    msf = panther.map { hmm, msf -> msf }
    PREPARE_TREEGRAFTER(
        SEARCH_PANTHER.out,
        msf
    )

    RUN_TREEGRAFTER(
        PREPARE_TREEGRAFTER.out.fasta,
        msf
    )

    PARSE_PANTHER(
        PREPARE_TREEGRAFTER.out.json.join(
            RUN_TREEGRAFTER.out, 
            by: [0, 1],
            failOnDuplicate: true,
            failOnMismatch: true
        )
    )
        
    emit:
    PARSE_PANTHER.out
        .map { meta, meta2, json -> tuple (meta, json) }  // [ meta, json ]
}