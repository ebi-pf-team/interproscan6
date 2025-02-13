process POPULATE_DATABASE {
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
    """
}

process INSERT_ORFS {
    // add protein seqs translated from nt seqs to the database
    input:
    val dbPath
    val fasta

    exec:
    println "TODO"
}

process BUILD_BATCHES {
    label 'local', 'ips6_container'

    input:
    val dbPath
    val batchSize

    output:
    tuple each tuple(val(meta), path(fasta))

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
