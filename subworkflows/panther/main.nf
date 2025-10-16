include { RUN_HMMER as SEARCH_PANTHER                         } from  "../../modules/hmmer"
include { PREPARE_TREEGRAFTER; RUN_TREEGRAFTER; PARSE_PANTHER } from  "../../modules/panther"

workflow PANTHER {
    take:
    ch_seqs
    dir
    hmm
    msf
    batch_size
    
    main:
    results = Channel.empty()

    SEARCH_PANTHER(
        ch_seqs,
        dir,
        hmm,
        "-Z 65000000 -E 0.001"
    )

    PREPARE_TREEGRAFTER(
        SEARCH_PANTHER.out,
        dir,
        msf
    )

    PREPARE_TREEGRAFTER.out.fasta
        .map { meta, sequenceIds, familyIds, fastas ->
            // Collate into batches of 500
            def batches = (0..<sequenceIds.size()).collect { i ->
                [sequenceIds[i], familyIds[i], fastas[i]]
            }.collate(500)
            // Return a list of tuples: [meta, [seqIds], [famIds], [fastas]]
            batches.withIndex().collect { batch, index ->
                def (seqs, fams, fas) = batch.transpose()
                tuple(meta, index, seqs, fams, fas)
            }
        }
        .flatten()
        .set { batched_fasta }

    batched_fasta.view()

    RUN_TREEGRAFTER(
        batched_fasta,
        dir,
        msf
    )

    results = results.mix(PARSE_PANTHER(
        PREPARE_TREEGRAFTER.out.json.join(RUN_TREEGRAFTER.out)
    ))
        
    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}