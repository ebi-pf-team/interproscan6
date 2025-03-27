process LOAD_SEQUENCES {
    // Populate a local sqlite3 database with sequences from the pipeline's input FASTA file.
    label         'local'
    errorStrategy 'terminate'

    input:
    val fasta
    val nucleic

    output:
    path "sequences.db"

    exec:
    def outputFilePath = task.workDir.resolve("sequences.db")
    SeqDB db = new SeqDB(outputFilePath.toString())
    db.loadFastaFile(fasta.toString(), nucleic)
    db.close()
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
