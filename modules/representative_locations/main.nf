import uk.ac.ebi.interpro.ProcessReprLocations

process REPRESENTATIVE_LOCATIONS_BULK {
    label    'mem_high', 'time_short', 'ips6_container'

    input:
    tuple val(meta), path(json)

    output:
    tuple val(meta), path("matches-with-repr-locations.json")

    script:
    """
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/bin/select-repr-locations.groovy ${json} matches-with-repr-locations.json
    """
}

process REPRESENTATIVE_LOCATIONS {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(json)

    output:
    tuple val(meta), path("matches-with-repr-locations.json")

    exec:
    def output = task.workDir.resolve('matches-with-repr-locations.json')
    ProcessReprLocations.run(json, output)
}
