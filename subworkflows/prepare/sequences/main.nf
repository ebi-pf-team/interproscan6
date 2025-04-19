include { VALIDATE_FASTA;
          LOAD_SEQUENCES;
          LOAD_ORFS;
          SPLIT_FASTA     } from "../../../modules/prepare_sequences"
include { ESL_TRANSLATE   } from "../../../modules/esl_translate"

workflow PREPARE_SEQUENCES {
    take:
    fasta_file   // path to the input FASTA file
    apps         // list of applications to be run

    main:
    def (validated_fasta, invalidCharCheckFile) = VALIDATE_FASTA(fasta_file, params.nucleic, apps, params.appsConfig)
    // Wait for process to complete and check the contents of the invalid_char_check file
    invalidCharCheckFile.subscribe { file ->
        def content = file.text
        if (content) {
            log.error "Invalid characters found in the input FASTA file:\n$content"
            exit 1
        }
    }

    if (params.nucleic) {
        // Store the input seqs in the internal ips6 seq db
        LOAD_SEQUENCES(validated_fasta, params.nucleic)

        // Chunk input file in smaller files for translation
        Channel.fromPath(params.input)
            .splitFasta( by: params.batchSize, file: true )
            .set { ch_fasta }

        /* Translate DNA/RNA sequences to protein sequences. Only proceed once completed
        ensuring SPLIT_FASTA only runs once LOAD_ORFS is completed */
        ch_translated = ESL_TRANSLATE(ch_fasta).collect()

        // Store sequences in the sequence database
        seq_db_path = LOAD_ORFS(ch_translated, LOAD_SEQUENCES.out)
    } else {
        // Store the input seqs in the internal ips6 seq db
        seq_db_path = LOAD_SEQUENCES(validated_fasta, params.nucleic)
    }
    // Build batches of unique protein seqs for the analysis
    SPLIT_FASTA(seq_db_path, params.batchSize, params.nucleic)

    fastaList = SPLIT_FASTA.out.collect()
    // Convert a list (or single file path) to a list of tuples containing indexed fasta file paths
    ch_seqs = fastaList
        .map { fastaList -> fastaList.indexed() } // creates a map-like object
        .flatMap()
        .map { entry -> [entry.key, entry.value] } // Convert to tuple [index, fasta]

    emit:
    ch_seqs      // a list of tuples: [index, fasta]
    seq_db_path  // str repr of path to the local sequence database
}
