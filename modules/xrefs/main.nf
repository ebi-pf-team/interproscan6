import java.nio.file.Path
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.ProcessXrefs

process XREFS_BULK {
    label    'mem_high', 'time_veryshort', 'ips6_container'

    input:
    tuple val(meta), path(json)
    path interpro_dir
    path panther_paint_dir
    val add_goterms
    val add_pathways

    output:
    tuple val(meta), path('matches-with-xrefs.json')

    script:
    """
    groovy -cp "/opt/interproscan6/lib:." /opt/interproscan6/bin/add-xrefs.groovy \
        ${json} \
        ${interpro_dir.name == 'DUMMY' ? '-' : interpro_dir} \
        ${panther_paint_dir.name == 'DUMMY2' ? '-' : panther_paint_dir} \
        ${add_goterms ? 'true' : 'false'} \
        ${add_pathways ? 'true' : 'false'} \
        matches-with-xrefs.json
    """
}

process XREFS {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(json)
    val interpro_dir
    val panther_paint_dir
    val add_goterms
    val add_pathways

    output:
    tuple val(meta), path('matches-with-xrefs.json')

    exec:
    def output = task.workDir.resolve('matches-with-xrefs.json')
    ProcessXrefs.run(json, 
                     interpro_dir.name == "DUMMY" ? null : interpro_dir, 
                     panther_paint_dir.name == "DUMMY2" ? null : panther_paint_dir,
                     add_goterms, add_pathways, output)
}
