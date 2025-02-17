process POPULATE_SEQ_DATABASE {
    // Populate a local sqlite3 database with sequences from the pipeline's input FASTA file.
    label 'local', 'ips6_container'

    input:
    val fasta
    val nucleic

    output:
    path "ips6.seq.db"

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py \
        ips6.seq.db \
        populate_sequences \
        --fasta $fasta \
         ${nucleic ? '--nucleic' : ''}
    chmod 777 ips6.seq.db
    """
}

process UPDATE_ORFS {
    // add protein seqs translated from ORFS in the nt seqs to the database
    maxForks 1  // Ensure that only one instance runs at a time to avoid concurrent writing to the db
    label 'local', 'ips6_container'

    input:
    val translatedFasta  // one FASTA per ESL_TRANSLATE batch
    val dbPath

    output:
    val dbPath

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py \
        $dbPath \
        update_orfs \
        --fasta $translatedFasta
    """
}

process BUILD_BATCHES {
    // Build the FASTA file batches of unique protein sequences for the sequence analysis
    label 'local', 'ips6_container'

    input:
    val dbPath
    val batchSize

    output:
    path("*.fasta")

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py \
        $dbPath \
        build_batches \
        --batch_size $batchSize
    """
}

process INDEX_FASTA_FILES {
    input:
    val fastaFiles

    output:
    val indexedFiles

    exec:
    indexedFiles = []
    int index = 1
    // handle when a single fasta file path is provided
    def fastaList = fastaFiles instanceof List ? fastaFiles : [fastaFiles]
    for (fasta: fastaList) {
        indexedFiles << [index, fasta]
        index += 1
    }
}
