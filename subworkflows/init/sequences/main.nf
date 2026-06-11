include { VALIDATE_FASTA;
          LOAD_SEQUENCES;
          LOAD_ORFS;
          SPLIT_FASTA     } from "../../../modules/prepare_sequences"
include { ESL_TRANSLATE   } from "../../../modules/esl_translate"

workflow INIT_SEQUENCES {
    take:
    fasta       // channel of the input FASTA file
    is_nucleic  // true if the sequences are nucleic, false if protein
    batch_size  // size of the batches to split the sequences into

    main:
    (ch_fasta, error) = VALIDATE_FASTA(fasta, is_nucleic)
    // Wait for process to complete so its output channels become available
    error.subscribe { seq_id -> 
        if (seq_id != null) {
            log.error "Invalid character(s) found in the input FASTA file."
            exit 1
        }
    }

    if (is_nucleic) {
        // Store the input seqs in the internal ips6 seq db
        LOAD_SEQUENCES(ch_fasta, is_nucleic)

        // Translate DNA/RNA sequences to protein sequences
        ESL_TRANSLATE(ch_fasta)

        // Store sequences in the sequence database
        seq_db = LOAD_ORFS(ESL_TRANSLATE.out, LOAD_SEQUENCES.out)
    } else {
        // Store the input seqs in the internal ips6 seq db
        seq_db = LOAD_SEQUENCES(ch_fasta, is_nucleic)
    }

    // Build batches of unique protein seqs for the analysis
    SPLIT_FASTA(seq_db, batch_size, is_nucleic)

    ch_split_fasta = SPLIT_FASTA.out.collect()

    // Convert a list (or single file path) to a list of tuples containing indexed fasta file paths
    ch_meta_fasta = ch_split_fasta
        .map { fasta_list -> fasta_list.indexed() } // creates a map-like object
        .flatMap()
        .map { entry -> [entry.key, entry.value] }  // Convert to tuple [index, fasta]

    emit:
    fasta = ch_meta_fasta   // a list of tuples: [index, fasta]
    database = seq_db       // path to the local sequence database
}
