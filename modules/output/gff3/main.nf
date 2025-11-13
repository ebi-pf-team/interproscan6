import uk.ac.ebi.interpro.ProcessOutputGFF3


process WRITE_GFF3 {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version

    exec:
    ProcessOutputGFF3.run(matches_files.collect { it.toString() }, 
                          seq_db_file.toString(), 
                          nucleic, 
                          interproscan_version, 
                          output_file.toString())
}
