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
    python3 $projectDir/interproscan/scripts/database.py ips6.seq.db $fasta populate_sequences $nucleic
    chmod 777 ips6.seq.db
    """
}

process UPDATE_ORFS {
    // add protein seqs translated from nt seqs to the database
    maxForks 1  // Ensure that only one instance runs at a time to avoid concurrent writing to the db
    label 'local', 'ips6_container'

    input:
    val translated_fasta  // one FASTA per ESL_TRANSLATE batch
    val db_path

    output:
    val ""  // to pass linting

    script:
    """
    python3 $projectDir/interproscan/scripts/database.py $db_path $translated_fasta update_orfs true > debug
    """
}

process BUILD_BATCHES {
    label 'local', 'ips6_container'

    input:
    val dbPath
    val batchSize

    output:
    val ""

    exec:
    // SequenceDatabase conn = new SequenceDatabase(dbPath)
    def offset = 0
    def meta = 1
    def files = []  // list to store results
//
//     while (true) {
//         def fasta = file("sequence.$meta.fasta")
//         fasta.text = ""
//
//         String query = """
//             SELECT protein_md5, sequence FROM PROTEIN_SEQUENCE
//             LIMIT $batchSize OFFSET $offset
//         """.toString()
//
//         def rows = conn.rows(query)
//         if (rows.isEmpty()) {
//             break
//         }
//         rows.each { row ->
//             fasta.append("> ${row.protein_md5}\n${row.sequence.replaceAll(/(.{1,60})/, '$1\n')}\n")
//         }
//         files << tuple(meta, fasta) // store the batch
//         offset += batchSize
//         meta += 1
//     }
//
//     return files
}
