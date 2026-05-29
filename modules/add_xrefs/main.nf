process ADD_XREFS {
    label     'mem_medium'
    label     'time_veryshort'
    container 'interpro/groovy:4.0.27-1'

    input:
    tuple val(meta), path(json, arity: '1..*', name: 'inputs/*')
    path interpro_dir
    path panther_paint_dir
    val add_goterms
    val add_pathways

    output:
    tuple val(meta), path('outputs/*.json', arity: '1..*')

    script:
    """
    add-xrefs.groovy \
        /opt/interproscan6/lib \
        inputs \
        ${interpro_dir.name == 'DUMMY' ? '-' : interpro_dir} \
        ${panther_paint_dir.name == 'DUMMY2' ? '-' : panther_paint_dir} \
        ${add_goterms ? 'true' : 'false'} \
        ${add_pathways ? 'true' : 'false'} \
        outputs

    chmod -R 777 outputs
    """
}
