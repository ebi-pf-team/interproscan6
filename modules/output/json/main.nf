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

    output:
    val output_file

    exec:
    ProcessOutputJSON.run(matches_files, 
                          seq_db_file, 
                          db_releases,
                          nucleic, 
                          interproscan_version,
                          jsonlines,
                          output_file)
}

process WRITE_JSON_BULK {
    label    'mem_veryhigh', 'time_long', 'ips6_container'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    path output_file
    path seq_db_file
    val nucleic
    val interproscan_version
    val db_releases
    val jsonlines

    output:
    val output_file

    script:
    def to_serialize = [:]
    db_releases.each { key, value ->
        to_serialize[key] = [
            version: value.version,
            dirpath: value.dirpath?.toString() 
        ]
    }
    def json = JsonOutput.toJson(to_serialize)
    """
    echo '${json}' > metadata.json
    groovy -cp "/opt/interproscan6/lib:/opt/interproscan6/lib/*:." /opt/interproscan6/bin/write-output.groovy \
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