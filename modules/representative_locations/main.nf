process REPRESENTATIVE_LOCATIONS {
    label     'mem_medium'
    label     'time_short'
    container 'interpro/groovy:4.0.27-1'

    input:
    tuple val(meta), path(json, arity: '1..*', name: 'inputs/*')

    output:
    tuple val(meta), path('outputs/*json', arity: '1..*')

    script:
    """
    select-repr-locations.groovy \
        /opt/interproscan6/lib \
        inputs \
        outputs

    chmod -R 777 outputs
    """
}
