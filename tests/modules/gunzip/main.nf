// codenarc-disable ProcessWithoutExecEnvironmentRule
process GUNZIP {
    cpus   1
    memory 1.GB

    input:
    path file_gz

    output:
    path "${file_gz.baseName}"

    script:
    """
    gunzip -c ${file_gz} > ${file_gz.baseName}
    """
}