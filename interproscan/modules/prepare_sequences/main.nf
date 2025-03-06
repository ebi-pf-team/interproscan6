process POPULATE_SEQ_DATABASE {
    // Populate a local sqlite3 database with sequences from the pipeline's input FASTA file.
    label         'local', 'ips6_container'
    errorStrategy 'terminate'

    input:
    file fasta
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
    label         'local', 'ips6_container'
    errorStrategy 'terminate'

    input:
    val translatedFastas  // could be one or multiple paths
    val dbPath

    output:
    val dbPath // ensure BUILD_BATCHES runs after UPDATE_ORFS

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py \
        $dbPath \
        update_orfs \
        --fasta "$translatedFastas"
    """
}

process BUILD_BATCHES {
    // Build the FASTA file batches of unique protein sequences for the sequence analysis
    label         'local', 'ips6_container'
    errorStrategy 'terminate'

    input:
    val dbPath
    val batchSize
    val nucleic

    output:
    path "*.fasta"

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py \
        $dbPath \
        build_batches \
        --batch_size $batchSize \
        ${nucleic ? '--nucleic' : ''}
    """
}

process INDEX_BATCHES {
    label         'local', 'ips6_container'
    errorStrategy 'terminate'

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
