import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessOutputJSON

process WRITE_JSON_LOCAL {
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

process WRITE_JSON {
    label    'mem_veryhigh', 'time_long', 'groovy_container'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    val output_file
    path seq_db_file
    val nucleic
    val interproscan_version
    val db_releases
    val jsonlines

    script:
    def json = JsonOutput.toJson(db_releases)
    """
    echo '${json}' > metadata.json
    groovy -cp "${projectDir}/lib:${projectDir}/lib/*:." ${projectDir}/bin/interpro/write-output.groovy \
        ${jsonlines ? 'jsonl' : 'json'} \
        . \
        metadata.json \
        ${seq_db_file} \
        ${nucleic ? 'true' : 'false'} \
        ${interproscan_version} \
        ${output_file}
    chmod 666 ${output_file}
    """
}