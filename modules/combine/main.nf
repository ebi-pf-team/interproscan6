import uk.ac.ebi.interpro.ProcessCombine

process COMBINE_MATCHES {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path('combined.json')

    exec:
    ProcessCombine.run(members_matches, task.workDir.resolve('combined.json'))
}

process COMBINE_MATCHES_BULK {
    label    'mem_high', 'time_veryshort', 'ips6_container'

    input:
    tuple val(meta), path(jsons, arity: '1..*', name: '?/*')

    output:
    tuple val(meta), path('combined.json')

    script:
    """
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/bin/combine.groovy . combined.json
    """
}