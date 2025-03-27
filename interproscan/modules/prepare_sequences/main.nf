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
    db.loadFastaFile(fasta.toString(), nucleic, false)
    db.close()
}

process LOAD_ORFS {
    // add protein seqs translated from ORFS in the nt seqs to the database
    label         'local'
    errorStrategy 'terminate'

    input:
    val translatedFastas  // could be one or multiple paths
    val dbPath

    output:
    val dbPath // ensure BUILD_BATCHES runs after LOAD_ORFS

    exec:
    SeqDB db = new SeqDB(dbPath.toString())
    translatedFastas.each {
        db.loadFastaFile(it.toString(), false, true)
    }
    db.close()
}

process SPLIT_FASTA {
    // Build the FASTA file batches of unique protein sequences for the sequence analysis
    label         'local'
    errorStrategy 'terminate'

    input:
    val dbPath
    val batchSize

    output:
    path "*.fasta"

    exec:
    String prefix = task.workDir.resolve("input").toString()
    SeqDB db = new SeqDB(dbPath.toString())
    db.splitFasta(prefix, batchSize)
    db.close()
}
