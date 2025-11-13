import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessOutputJSON

process WRITE_JSON {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matches_files  // {query prot seq md5: {model acc: match}}
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val db_releases
    val jsonlines

    exec:
    ProcessOutputJSON.run(matches_files.collect { it.toString() }, 
                          seq_db_file.toString(), 
                          db_releases,
                          nucleic, 
                          interproscan_version,
                          jsonlines,
                          output_file.toString())
}
