import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessXrefs

process XREFS {
    label    'mem_high', 'time_veryshort', 'groovy_container'

    input:
    tuple val(meta), path(json_combined)
    val db_releases
    val add_goterms
    val add_pathways
    val panther_paint_dir

    output:
    tuple val(meta), path('matches2xrefs.json')

    script:
    def json = JsonOutput.toJson(db_releases)
    """
    echo '${json}' > metadata.json
    groovy -cp "${projectDir}/lib:." ${projectDir}/bin/interpro/add-xrefs.groovy \
        ${json_combined} \
        metadata.json \
        ${add_goterms ? 'true' : 'false'} \
        ${add_pathways ? 'true' : 'false'} \
        ${panther_paint_dir} \
        matches2xrefs.json
    """
}

process XREFS_LOCAL {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(matches_path)
    val db_releases
    val add_goterms
    val add_pathways
    val panther_paint_dir

    output:
    tuple val(meta), path('matches2xrefs.json')

    exec:
    String outputFilePath = task.workDir.resolve('matches2xrefs.json')
    ProcessXrefs.run(matches_path.toString(), db_releases, add_goterms, add_pathways, panther_paint_dir, outputFilePath)
}

