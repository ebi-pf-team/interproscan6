import uk.ac.ebi.interpro.ProcessOutputTSV

process WRITE_TSV_LOCAL {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matches_files
    val output_file
    val seq_db_file
    val nucleic

    output:
    val output_file

    exec:
    ProcessOutputTSV.run(
        matches_files.collect { it.toString() },
        seq_db_file.toString(),
        nucleic,
        output_file.toString()
    )
}

process WRITE_TSV {
    label    'mem_veryhigh', 'time_long', 'groovy_container'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    val output_file
    path seq_db_file
    val nucleic

    output:
    val output_file

    script:
    """
    groovy -cp "${projectDir}/lib:${projectDir}/lib/*:." ${projectDir}/bin/interpro/write-output.groovy \
        tsv \
        . \
        - \
        ${seq_db_file} \
        ${nucleic ? 'true' : 'false'} \
        - \
        ${output_file}
    chmod 666 ${output_file}
    """
}
