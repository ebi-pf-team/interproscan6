import uk.ac.ebi.interpro.ProcessOutputTSV

process WRITE_TSV {
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
        matches_files,
        seq_db_file,
        nucleic,
        output_file
    )
}

process WRITE_TSV_BULK {
    label     'mem_veryhigh', 'time_long'
    container 'interpro/groovy:4.0.27-1'

    input:
    path(input_files, arity: '1..*', name: '?/*')
    path output_file
    path seq_db_file
    val nucleic

    output:
    path output_file

    script:
    """
    groovy -cp "/opt/interproscan6/lib:/opt/interproscan6/lib/*:." /opt/interproscan6/write-output.groovy \
        tsv \
        . \
        ${seq_db_file} \
        ${nucleic ? 'true' : 'false'} \
        - \
        - \
        ${output_file}
    chmod 666 ${output_file}
    """
}
