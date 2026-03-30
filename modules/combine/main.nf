import uk.ac.ebi.interpro.ProcessCombine

process COMBINE_MATCHES_LOCAL {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path('combined.json')

    exec:
    ProcessCombine.run(members_matches, task.workDir.resolve('combined.json'))
}

process COMBINE_MATCHES {
    label    'mem_high', 'time_veryshort', 'groovy_container'

    input:
    tuple val(meta), path(jsons, arity: '1..*', name: '?/*')

    output:
    tuple val(meta), path('combined.json')

    script:
    """
    groovy -cp "${projectDir}/lib:." ${projectDir}/bin/interpro/combine.groovy . combined.json
    """
}