import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.SeqDB

process VALIDATE_FASTA {
    // check the formating of the intput FASTA, i.e. look for illegal characters
    label         'mem_low', 'time_short'
    executor      'local'
    errorStrategy 'terminate'

    input:
    val fasta
    val is_nucleic

    output:
    val fasta
    val seq_id

    exec:
    seq_id = FastaFile.validate(fasta.toString(), is_nucleic)
}

process LOAD_SEQUENCES {
    // Populate a native sqlite3 database with sequences from the pipeline's input FASTA file.
    label         'mem_low', 'time_short'
    executor      'local'
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
    label         'mem_low', 'time_short'
    executor      'local'
    errorStrategy 'terminate'

    input:
    val translated_fasta
    val db_path

    output:
    val db_path // ensure BUILD_BATCHES runs after LOAD_ORFS

    exec:
    SeqDB db = new SeqDB(db_path.toString())
    db.loadFastaFile(translated_fasta.toString(), false, true)
    db.close()
}

process SPLIT_FASTA {
    // Build the FASTA file batches of unique protein sequences for the sequence analysis
    label         'mem_low', 'time_short'
    executor      'local'
    errorStrategy 'terminate'

    input:
    val db_path
    val batch_size
    val nucleic

    output:
    path "*.fasta"

    exec:
    String prefix = task.workDir.resolve("input").toString()
    SeqDB db = new SeqDB(db_path.toString())
    db.splitFasta(prefix, batch_size, nucleic)
    db.close()
}