include { RUN_HMMER as SEARCH_GENE3D                                  } from  "../../../modules/hmmer"
include { RESOLVE_GENE3D; ASSIGN_CATH; PARSE_CATHGENE3D               } from  "../../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; RESOLVE_FUNFAM; PARSE_FUNFAM } from  "../../../modules/cath/funfam"

workflow CATH {
    take:
    ch_seqs
    applications
    gene3d_hmm
    gene3d_model2sfs
    gene3d_disc_regs
    funfam_dir
    funfam_chunksize

    main:
    results = Channel.empty()

    // Search Gene3D profiles
    SEARCH_GENE3D(
        ch_seqs,
        gene3d_hmm,
        "-Z 65245 -E 0.001"
    )
    ch_gene3d = SEARCH_GENE3D.out

    // Select best domain matches
    RESOLVE_GENE3D(ch_gene3d)

    // Assign CATH superfamily to matches
    ASSIGN_CATH(
        RESOLVE_GENE3D.out,
        gene3d_model2sfs,
        gene3d_disc_regs
    )

    // Join results and parse them
    ch_gene3d = ch_gene3d.join(ASSIGN_CATH.out)
    PARSE_CATHGENE3D(ch_gene3d)

    if (applications.contains("gene3d")) {
        results = results.mix(PARSE_CATHGENE3D.out)
    }

    if (applications.contains("cathfunfam")) {
        // Find unique CATH superfamilies with at least one hit
        PREPARE_FUNFAM(
            PARSE_CATHGENE3D.out,
            funfam_dir
        )

        /*
            Join input fasta file with superfamilies.
            We split in smaller chunks to parallelize searching FunFam profiles
        */
        ch_seqs
            .join(PREPARE_FUNFAM.out)
            .flatMap { id, file, supfams ->
                supfams
                    .collate(funfam_chunksize)
                    .collect { chunk -> [id, file, chunk] }
            }
            .set { ch_funfams }

        // Search FunFam profiles
        SEARCH_FUNFAM(
            ch_funfams,
            funfam_dir
        )

        // Select best domain matches
        RESOLVE_FUNFAM(SEARCH_FUNFAM.out)

        // Join results and parse them
        ch_cathfunfam = SEARCH_FUNFAM.out.join(RESOLVE_FUNFAM.out)
        PARSE_FUNFAM(ch_cathfunfam)
        results = results.mix(PARSE_FUNFAM.out)
    }

    emit:
    results
}