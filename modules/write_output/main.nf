process WRITE_GFF3 {
    label    'mem_medium'
    label    'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val interpro_version

    output:
    val output_file

    exec:
    def localDbPath = null
    try {
        localDbPath = uk.ac.ebi.interpro.LocalSeqDB.copyFrom(seq_db_file, params.localSeqDbDir)
        uk.ac.ebi.interpro.ProcessOutputGFF3.run(
            matches_files,
            localDbPath,
            nucleic,
            interproscan_version,
            interpro_version,
            output_file)
    } finally {
        uk.ac.ebi.interpro.LocalSeqDB.cleanup(localDbPath)
    }
}

process WRITE_JSON {
    label    'mem_medium'
    label    'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val interpro_version
    val jsonlines

    output:
    val output_file

    exec:
    def localDbPath = null
    try {
        localDbPath = uk.ac.ebi.interpro.LocalSeqDB.copyFrom(seq_db_file, params.localSeqDbDir)
        uk.ac.ebi.interpro.ProcessOutputJSON.run(
            matches_files,
            localDbPath,
            nucleic,
            interproscan_version,
            interpro_version,
            jsonlines,
            output_file)
    } finally {
        uk.ac.ebi.interpro.LocalSeqDB.cleanup(localDbPath)
    }
}

process WRITE_TSV {
    label    'mem_medium'
    label    'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic

    output:
    val output_file

    exec:
    def localDbPath = null
    try {
        localDbPath = uk.ac.ebi.interpro.LocalSeqDB.copyFrom(seq_db_file, params.localSeqDbDir)
        uk.ac.ebi.interpro.ProcessOutputTSV.run(
            matches_files,
            localDbPath,
            nucleic,
            output_file)
    } finally {
        uk.ac.ebi.interpro.LocalSeqDB.cleanup(localDbPath)
    }
}

process WRITE_XML {
    label    'mem_medium'
    label    'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val interpro_version

    output:
    val output_file

    exec:
    def localDbPath = null
    try {
        localDbPath = uk.ac.ebi.interpro.LocalSeqDB.copyFrom(seq_db_file, params.localSeqDbDir)
        uk.ac.ebi.interpro.ProcessOutputXML.run(
            matches_files,
            localDbPath,
            nucleic,
            interproscan_version,
            interpro_version,
            output_file)
    } finally {
        uk.ac.ebi.interpro.LocalSeqDB.cleanup(localDbPath)
    }
}   
