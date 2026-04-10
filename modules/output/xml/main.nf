import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessOutputXML

process WRITE_XML {
    label    'mem_low', 'time_long'
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
    ProcessOutputXML.run(matches_files, 
                         seq_db_file, 
                         nucleic, 
                         interproscan_version,
                         interpro_version,
                         output_file)
}

process WRITE_XML_BULK {
    label      'mem_veryhigh', 'time_long'
    container 'interpro/groovy:4.0.27-1'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    path output_file
    path seq_db_file
    val nucleic
    val interproscan_version
    val interpro_version

    output:
    path output_file

    script:
    """
    groovy -cp "/opt/interproscan6/lib:/opt/interproscan6/lib/*:." /opt/interproscan6/write-output.groovy \
        xml \
        . \
        ${seq_db_file} \
        ${nucleic ? 'true' : 'false'} \
        ${interproscan_version} \
        ${interpro_version ?: '-'} \
        ${output_file}
    chmod 666 ${output_file}
    """
}