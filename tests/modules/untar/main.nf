// codenarc-disable ProcessWithoutExecEnvironmentRule
process UNTAR {
    cpus   1
    memory 1.GB

    input:
    path archive

    output:
    path "${directory}"

    script:
    directory = archive.simpleName
    """
    mkdir ${directory}

    tar -C ${directory} -xavf ${archive} 
    """
}