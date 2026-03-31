import uk.ac.ebi.interpro.ProcessReprLocations

process REPRESENTATIVE_LOCATIONS_BULK {
    label    'mem_high', 'time_short', 'ips6_container'

    input:
    tuple val(meta), path(matches_path)

    output:
    tuple val(meta), path("matches_with_representative.json")

    script:
    """
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/bin/select-repr-locations.groovy ${matches_path.toString()} matches_with_representative.json
    """
}

process REPRESENTATIVE_LOCATIONS {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(matches_path)

    output:
    tuple val(meta), path("matches_with_representative.json")

    exec:
    Path outputFilePath = task.workDir.resolve('matches_with_representative.json')
    ProcessReprLocations.run(matches_path, outputFilePath)
}
