import uk.ac.ebi.interpro.ProcessReprLocations

process REPRESENTATIVE_LOCATIONS {
    label    'mem_high', 'time_short', 'groovy_container'

    input:
    tuple val(meta), val(matches_path)

    output:
    tuple val(meta), path("matches_with_representative.json")

    script:
    """
    groovy -cp "${projectDir}/lib:." ${projectDir}/bin/interpro/select-repr-locations.groovy ${matches_path.toString()} matches_with_representative.json
    """
}

process REPRESENTATIVE_LOCATIONS_LOCAL {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(matches_path)

    output:
    tuple val(meta), path("matches_with_representative.json")

    exec:
    String outputFilePath = task.workDir.resolve('matches_with_representative.json')
    ProcessReprLocations.run(matches_path.toString(), outputFilePath)
}
