import uk.ac.ebi.interpro.ProcessOutputTSV

process WRITE_TSV {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic

    exec:
    ProcessOutputTSV.run(
        matches_files.collect { it.toString() },
        seq_db_file.toString(),
        nucleic,
        output_file.toString()
    )
}

