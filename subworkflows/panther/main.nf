include { SEARCH_PANTHER; PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "${moduleDir}/../../modules/panther"

workflow PANTHER {
    take:
    ch_seqs    // channel of tuples (meta, fasta file)
    dir        // str repr of the data directory path
    hmm        // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    msf        // str repr of the path to the MSF file in the data dir -> datadir/msfFile
    batch_size // int, number of sequences per sub batch for searching

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
        PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out, by: [0, 1])
    ))
        
    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}