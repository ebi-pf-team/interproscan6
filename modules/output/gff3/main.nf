import uk.ac.ebi.interpro.ProcessOutputGFF3


process WRITE_GFF3_LOCAL {
    label    'mem_low', 'time_long'
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
    ProcessOutputGFF3.run(matches_files.collect { it.toString() }, 
        seq_db_file.toString(), 
        nucleic, 
        interproscan_version, 
        output_file.toString())
}

process WRITE_GFF3 {
    label    'mem_veryhigh', 'time_long', 'groovy_container'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    val output_file
    path seq_db_file
    val nucleic
    val interproscan_version

    output:
    val output_file

    script:
    """
    groovy -cp "${projectDir}/lib:${projectDir}/lib/*:." ${projectDir}/bin/interpro/write-output.groovy \
        gff3 \
        . \
        - \
        ${seq_db_file} \
        ${nucleic ? 'true' : 'false'} \
        ${interproscan_version} \
        ${output_file}
    chmod 666 ${output_file}
    """
}
