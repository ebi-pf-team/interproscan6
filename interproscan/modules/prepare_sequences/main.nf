import groovy.json.JsonOutput
import groovy.sql.Sql

process POPULATE_DATABASE {
    label 'local'

    input:
    val fasta
    val dbPath
    val nucleic

    output:
    val "" // to pass linting, because all processes have to have an output

    exec:
    SequenceDatabase conn = SequenceDatabase(dbPath)
    def currentSequence = null
    new File(fasta).eachLine { line ->
        if (line.startswith(">")) {
            if (currentSequence) {
                if (nucleic) {
                    addNucleotideToDb(currentSequence, conn)
                } else {
                    addProteinToDb(currentSequence, conn)
                }
            }
            def header = line.substring(1).split(" ", 2).collect { it.trim() }
            currentSequence = [
                id          : header[0],
                description : header[1],
                sequence    : ""
            ]
        } else {
            currentSequence.sequence += line.trim()
        }
    }
    if (currentSequence) { // make sure to add the last seq in the file
        if (nucleic) {
            addNucleotideToDb(currentSequence, conn)
        } else {
            addProteinToDb(currentSequence, conn)
        }
    }
    conn.close()
}

def addNucleotideToDb(Map currentSequence, SequenceDatabase conn) {
    currentSequence["md5"] = FastaSequence.getMD5(currentSequence.sequence)
    conn.insertNucleotideSequence(currentSequence.md5, currentSequence.sequence)
    conn.insertNucleotide(currentSequence.md5, currentSequence.id, currentSequence.description)
}

def addProteinToDb(Map currentSequence, SequenceDatabase conn) {
    currentSequence["md5"] = FastaSequence.getMD5(currentSequence.sequence)
    conn.insertProteinSequence(currentSequence.md5, currentSequence.sequence)
    conn.insertProtein(currentSequence.md5, currentSequence.id, currentSequence.description)
}

process UPDATE_DATABASE {
    // add protein seqs translated from nt seqs to the database
    input:
    val dbPath
    val fasta

    output:
    val "" // to pass linting, because all processes have to have an output

    exec:
    println "TODO"
}

process BUILD_BATCHES {
    label 'local'

    input:
    val dbPath
    val batchSize

    output:
    tuple each tuple(val(meta), path(fasta))

    exec:
    SequenceDatabase conn = SequenceDatabase(dbPath)
    def offset = 0
    int meta = 1
    List<Tuple> files = []  // list to store results

    while (true) {
        def fasta = file("sequence.$meta.fasta")
        fasta.text = ""

        String query = """
            SELECT protein_md5, sequence FROM PROTEIN_SEQUENCE
            LIMIT $batchSize OFFSET $offset
        """.toString()

        def rows = conn.rows(query)
        if (rows.isEmpty()) {
            break
        }
        rows.each { row ->
            fasta.append("> ${row.protein_md5}\n${row.sequence.replaceAll(/(.{1,60})/, '$1\n')}\n")
        }
        files << tuple(meta, fasta) // store the batch
        offset += batchSize
        meta += 1
    }

    emit: files
}
