include { SEARCH_GENE3D; PARSE_CATHGENE3D } from  "${moduleDir}/../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; PARSE_FUNFAM  } from  "${moduleDir}/../../modules/cath/funfam"

workflow CATH {
    take:
    ch_seqs               // channel of tuples (index, fasta file)
    report_cathgene3d     // boolean to report Gene3D results
    cathgene3d_dir        // str repr of the data directory path for Gene3D
    cathgene3d_hmm        // str repr of the path to the HMM file in the data dir -> datadir/hmmFile
    cathgene3d_model2sfs  // str repr of the path to the model2sf file in the data dir -> datadir/model2sfsFile
    cathgene3d_disc_regs  // str repr of the path to the disordered regions file in the data dir -> datadir/discRegsFile
    report_cathfunfam     // boolean to report FunFam results
    cathfunfam_dir        // str repr of the data directory path for FunFam
    batch_size            // int, number of sequences per sub batch for searching

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
    
    // Search Gene3D profiles, selec
    /* Select Gene3D profiles,
       select the best domain matches,
       and assign CATH superfamily to matches */
    SEARCH_GENE3D(
        ch_split,
        cathgene3d_dir,
        cathgene3d_hmm,
        cathgene3d_model2sfs,
        cathgene3d_disc_regs
    )
       
    PARSE_CATHGENE3D(SEARCH_GENE3D.out)

    if (report_cathgene3d) {
        results = results.mix(PARSE_CATHGENE3D.out)
    }

    if (report_cathfunfam) {
        // Find unique CATH superfamilies with at least one hit
        PREPARE_FUNFAM(
            PARSE_CATHGENE3D.out,
            cathfunfam_dir
        )

        // Join input fasta file with superfamilies.
        ch_funfams = ch_split
            .join(PREPARE_FUNFAM.out, by: [0, 1])

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