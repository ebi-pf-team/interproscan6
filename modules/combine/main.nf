import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

import uk.ac.ebi.interpro.Process

process COMBINE_MATCHES_LOCAL {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path('combined.json')

    exec:
    String outputFilePath = task.workDir.resolve('combined.json')
    Process.combineMatches(members_matches.collect { it.toString() }, outputFilePath)
}

process COMBINE_MATCHES {
    label    'mem_high', 'time_medium', 'groovy_container'

    input:
    tuple val(meta), path(jsons, arity: '1..*', name: '?/*')

    output:
    tuple val(meta), path('combined.json')

    script:
    """
    groovy -cp "${projectDir}/lib:." ${projectDir}/bin/interpro/combine.groovy . combined.json
    """
}