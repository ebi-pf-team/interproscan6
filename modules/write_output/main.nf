process WRITE_GFF3 {
    label    'mem_medium', 'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version

    output:
    val output_file

    exec:
    uk.ac.ebi.interpro.ProcessOutputGFF3.run(
        matches_files, 
        seq_db_file, 
        nucleic, 
        interproscan_version, 
        output_file)
}

process WRITE_JSON {
    label    'mem_medium', 'time_long'
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
    uk.ac.ebi.interpro.ProcessOutputJSON.run(
        matches_files, 
        seq_db_file, 
        nucleic, 
        interproscan_version,
        interpro_version,
        jsonlines,
        output_file)
}

process WRITE_TSV {
    label    'mem_medium', 'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic

    output:
    val output_file

    exec:
    uk.ac.ebi.interpro.ProcessOutputTSV.run(
        matches_files,
        seq_db_file,
        nucleic,
        output_file)
}

process WRITE_XML {
    label    'mem_medium', 'time_long'
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
    uk.ac.ebi.interpro.ProcessOutputXML.run(
        matches_files, 
        seq_db_file, 
        nucleic, 
        interproscan_version,
        interpro_version,
        output_file)
}   
