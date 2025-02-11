import groovy.json.JsonOutput
import groovy.sql.Sql

process PREPARE_NUCLEIC_SEQUENCES {
    label 'local'

    input:
    tuple val(n_fasta), val(p_fasta)

    output:
    tuple val(task.index), val(p_fasta), path("sequences.json")

    exec:
    def sequences = FastaFile.parse(n_fasta.toString())
    def ntSequences = [:]
    sequences.each { seq ->
        ntSequences[seq.id] = seq
    }

    sequences = FastaFile.parse(p_fasta.toString())
    def output = [:]
    sequences.each { seq ->
        def match = seq.description =~ /source=(\S+)/
        assert match
        def source = match[0][1]
        def ntSeq = ntSequences[source]
        seq.translatedFrom = ntSeq
        if (!output.containsKey(seq.md5)) {
            output[seq.md5] = []
        }
        output[seq.md5] << seq
    }
    def outputPath = task.workDir.resolve("sequences.json")
    def json = JsonOutput.toJson(output)
    new File(outputPath.toString()).write(json)
}

process NUCLEOTIDE_SEQUENCES_INTO_DB {
    label 'local'

    input:
    val fasta
    val dbPath

    output:
    val "" // to pass linting, because all processes have to have an output

    exec:
    SequenceDatabase conn = SequenceDatabase(dbPath)
    def currentSequence = null
    new File(fasta).eachLine { line ->
        if (line.startswith(">")) {
            if (currentSequence) {
                addNucleotideToDb(currentSequence, conn)
            }
            def header = line.substring(1).split(" ", 2).collect { it}
        }
    }
}

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
