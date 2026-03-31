import java.nio.file.Path
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessXrefs

process XREFS_BULK {
    label    'mem_high', 'time_veryshort', 'ips6_container'

    input:
    tuple val(meta), path(json_combined)
    val db_releases
    val add_goterms
    val add_pathways
    val panther_paint_dir

    output:
    tuple val(meta), path('matches2xrefs.json')

    script:
    def to_serialize = [:]
    db_releases.each { key, value ->
        to_serialize[key] = [
            version: value.version,
            dirpath: value.dirpath?.toString() 
        ]
    }
    def json = JsonOutput.toJson(to_serialize)
    """
    echo '${json}' > metadata.json
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/bin/add-xrefs.groovy \
        ${json_combined} \
        metadata.json \
        ${add_goterms ? 'true' : 'false'} \
        ${add_pathways ? 'true' : 'false'} \
        ${panther_paint_dir} \
        matches2xrefs.json
    """
}

process XREFS {
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
    def output = task.workDir.resolve('matches2xrefs.json')

    /*
    nf-test doesn't seem to support casting nested objects to file(), so we keep them a strings and convert them here
    */
    db_releases.each { _, value ->
        if (value.dirpath && value.dirpath.getClass().equals(String)) {
            value.dirpath = file(value.dirpath)
        }
    }

    ProcessXrefs.run(matches_path, db_releases, add_goterms, add_pathways, panther_paint_dir, output)
}
