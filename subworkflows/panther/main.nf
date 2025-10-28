include { SEARCH_PANTHER; PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "../../modules/panther"

workflow PANTHER {
    take:
    ch_seqs
    dir
    hmm
    msf
    batch_size
    
    main:
    results = Channel.empty()

    ch_split = ch_seqs
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()

    SEARCH_PANTHER(
        ch_split,
        dir,
        hmm
    )

    PREPARE_TREEGRAFTER(
        SEARCH_PANTHER.out,
        dir,
        msf
    )

    RUN_TREEGRAFTER(
        PREPARE_TREEGRAFTER.out.fasta,
        dir,
        msf
    )

    results = results.mix(PARSE_PANTHER(
        PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out)
    ))
        
    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}