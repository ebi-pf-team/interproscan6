import java.io.File
import java.nio.file.*
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

import uk.ac.ebi.interpro.InterProScan

process DOWNLOAD {
    maxForks 1
    executor 'local'
    label    'ips6_container'

    input:
    tuple val(name), val(arcname), val(version), val(skip), val(dirpath)
    val iprscan_version
    val use_globus
    path outdir

    output:
    tuple val(name), val(dirpath)

    script:
    if (skip) {
        """
        """
    } else {
        def base_url = use_globus ? InterProScan.GLOBUS_URL : InterProScan.FTP_URL
        """
        curl -OJ ${base_url}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz
        curl -OJ ${base_url}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz.md5
        md5sum -c ${arcname}-${version}.tar.gz.md5 || { echo "Error: MD5 checksum failed" >&2; exit 1; }
        tar -C ${outdir} -zxf ${arcname}-${version}.tar.gz
        rm ${arcname}-${version}.tar.gz*
        chmod 777 -R ${outdir}/${arcname}
        """        
    }
    
}

process FIND_DATABASES {
    label    'mem_min', 'time_veryshort'
    executor 'local'

    input:
    tuple val(n), val(p)  // state dependency (str, str)
    val databases_json    // path      
    val applications      // list[str]
    val appl_dirs         // map
    val datadir           // path

    output:
    val ready,    emit: ready
    val missing,  emit: missing

    exec:
    def json = new JsonSlurper().parse(databases_json)
    def normalised_json = [:]
    json.each { key, value ->
        normalised_json[key.replaceAll(/[\s\-]+/, '').toLowerCase()] = value
    }

    /* 
        Ugly trick to make sure that metadata for both CATH-Gene3D are CATH-FunFam
        are available even if we need only one
    */
    def appls = applications.clone() as Set
    if (appls.contains("cathgene3d") || appls.contains("cathfunfam")) {
        appls.add("cathgene3d")
        appls.add("cathfunfam")
    }

    ready   = [] as Set
    missing = [] as Set
    seen    = [] as Set
    appls.each { appl_name ->
        if (appl_dirs.containsKey(appl_name)) {
            def normalised_name = appl_name.replaceAll(/[\s\-]+/, '').toLowerCase()
            assert normalised_json.containsKey(normalised_name)
            def appl_version = normalised_json[normalised_name]
            def appl_dir_parts = appl_dirs[normalised_name]
            def appl_dir = appl_dir_parts[0]
            def appl_subdir = appl_dir_parts[1]
            def appl_fulldir = datadir.resolve("${appl_dir}/${appl_version}/${appl_subdir}")
            if (appl_fulldir.isDirectory()) {
                ready.add( [ normalised_name, appl_fulldir ])
            } else {
                missing.add( [ 
                    normalised_name, 
                    appl_dir, 
                    appl_version, 
                    seen.contains(appl_dir),
                    appl_fulldir 
                ] )
                seen.add( appl_dir )
            }
        }
    }
}

process COLLECT_DATABASES {
    label    'mem_min', 'time_veryshort'
    executor 'local'
    cache false  // Stops the esotericsoftware.kryo.serializers warning

    input:
    val list_databases

    output:
    val map_databases

    exec:
    map_databases = list_databases.collectEntries { entry ->
        def name = entry[0]
        def path = entry[1]
        [(name): path]
    }
}
