import uk.ac.ebi.interpro.ProcessCombine

process COMBINE_MATCHES {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(json)

    output:
    tuple val(meta), path('combined.json')

    exec:
    ProcessCombine.run(json, task.workDir.resolve('combined.json'))
}

process COMBINE_MATCHES_BULK {
    label     'mem_high', 'time_veryshort'
    container 'interpro/groovy:4.0.27-1'

    input:
    tuple val(meta), path(json, arity: '1..*', name: '?/*')

    output:
    tuple val(meta), path('combined.json')

    script:
    """
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/combine.groovy . combined.json
    """
}