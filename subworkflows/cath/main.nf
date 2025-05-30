include { SEARCH_GENE3D; RESOLVE_GENE3D; ASSIGN_CATH; PARSE_CATHGENE3D } from  "../../modules/cath/gene3d"
include { PREPARE_FUNFAM; SEARCH_FUNFAM; RESOLVE_FUNFAM; PARSE_FUNFAM  } from  "../../modules/cath/funfam"

workflow CATH {
    take:
    ch_seqs
    report_cathgene3d
    cathgene3d_dir
    cathgene3d_hmm
    cathgene3d_model2sfs
    cathgene3d_disc_regs
    report_cathfunfam
    cathfunfam_dir

    main:
    results = Channel.empty()

    ch_split = ch_seqs
        .map { meta, fasta ->
            fasta
                .splitFasta( by: 1000, file: true )
                .indexed()
                .collect { index, chunk -> [meta, index, chunk] }
        }
        .flatMap()
    
    // Search Gene3D profiles
    SEARCH_GENE3D(
        ch_split,
        cathgene3d_dir,
        cathgene3d_hmm
    )
    ch_gene3d = SEARCH_GENE3D.out

    // Select best domain matches
    RESOLVE_GENE3D(ch_gene3d)

    // Assign CATH superfamily to matches
    ASSIGN_CATH(
        RESOLVE_GENE3D.out,
        cathgene3d_dir,
        cathgene3d_model2sfs,
        cathgene3d_disc_regs
    )

    // Join results and parse them
    ch_cathgene3d = ch_gene3d
        .join(ASSIGN_CATH.out, by: [0, 1])
   
    PARSE_CATHGENE3D(ch_cathgene3d)

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

        // Search FunFam profiles
        SEARCH_FUNFAM(
            ch_funfams,
            cathfunfam_dir
        )

        // Select best domain matches
        RESOLVE_FUNFAM(SEARCH_FUNFAM.out)

        // Join results and parse them
        ch_cathfunfam = SEARCH_FUNFAM.out.join(RESOLVE_FUNFAM.out, by: [0, 1])
        PARSE_FUNFAM(ch_cathfunfam)
        results = results.mix(PARSE_FUNFAM.out)
    }

    emit:
    results.map { meta, meta2, json -> tuple (meta, json) }
}