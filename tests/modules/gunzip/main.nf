process GUNZIP {
    input:
    path file_gz

    output:
    path "${file_gz.simpleName}"

    script:
    """
    gunzip -c ${file_gz} > ${file_gz.simpleName}
    """
}