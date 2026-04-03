process UNTAR {
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