include { SEARCH_GENE3D; PARSE_CATHGENE3D } from  "../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; PARSE_FUNFAM  } from  "../../modules/cath/funfam"

workflow CATH {
    take:
    fasta                 // [meta, fasta]
    report_cathgene3d     // boolean to report Gene3D results
    cathgene3d_files      // [hmm, model2sfs, disc_regs]
    report_cathfunfam     // boolean to report FunFam results
    cathfunfam_dir        // path to the data directory path for FunFam
    batch_size            // int, number of sequences per sub batch for searching

    main:
    results = Channel.empty()

    ch_split = fasta
        .map { meta, fasta ->
            fasta
                .splitFasta( by: batch_size, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()
    
    /* Select Gene3D profiles,
       select the best domain matches,
       and assign CATH superfamily to matches */
    SEARCH_GENE3D(
        ch_split,
        cathgene3d_files
    )
       
    PARSE_CATHGENE3D(SEARCH_GENE3D.out)

    if (report_cathgene3d) {
        results = results.mix(PARSE_CATHGENE3D.out)
    }

    if (report_cathfunfam) {
        // Find unique CATH superfamilies with at least one hit
        PREPARE_FUNFAM(
            PARSE_CATHGENE3D.out
        )

        // Join input fasta file with superfamilies.
        ch_funfams = ch_split
            .join(
                PREPARE_FUNFAM.out, 
                by: [0, 1],
                failOnDuplicate: true,
                failOnMismatch: true
            )

        // Search FunFam profiles, and select the best domain matches
        SEARCH_FUNFAM(
            ch_funfams,
            cathfunfam_dir
        )

        // Parse results
        PARSE_FUNFAM(SEARCH_FUNFAM.out)
        results = results.mix(PARSE_FUNFAM.out)
    }

    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}