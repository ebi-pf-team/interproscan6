import java.io.File
import java.nio.file.*
import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process DOWNLOAD {
    maxForks 1
    executor 'local'
    label    'ips6_container'

    input:
    tuple val(name), val(arcname), val(version), val(skip), val(path)
    val iprscan_version
    path outdir

    output:
    tuple val(name), val(version), val(path)

    script:
    if (skip) {
        """
        """
    } else {
        """
        curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz
        curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${arcname}/${arcname}-${version}.tar.gz.md5
        md5sum -c ${arcname}-${version}.tar.gz.md5 || { echo "Error: MD5 checksum failed" >&2; exit 1; }
        tar -C ${outdir} -zxf ${arcname}-${version}.tar.gz
        rm ${arcname}-${version}.tar.gz*
        chmod 777 -R ${outdir}/${arcname}
        """        
    }
    
}

process FIND_MISSING_DATA {
    label    'tiny'
    executor 'local'

    input:
    tuple val(n), val(v), val(p)  // state dependency
    val json_database
    val apps_to_run
    val app_dirs
    val datadir

    output:
    val with_data,          emit: with_data
    val without_data,       emit: without_data

    exec:
    File file = new File(json_database)
    def json = new JsonSlurper().parse(file)
    def normalised_json = [:]
    json.each { key, value ->
        normalised_json[key.replaceAll(/[\s\-]+/, '').toLowerCase()] = value
    }

    /* 
        Ugly trick to make sure that metadata for both CATH-Gene3D are CATH-FunFam
        are available even if we need only one
    */
    if (apps_to_run.contains("cathgene3d") || apps_to_run.contains("cathfunfam")) {
        apps_to_run.add("cathgene3d")
        apps_to_run.add("cathfunfam")
    }

    with_data = [] as Set
    without_data = [] as Set
    to_download = [] as Set
    apps_to_run.each { db_name ->
        if (app_dirs.containsKey(db_name)) {
            def normalised_name = db_name.replaceAll(/[\s\-]+/, '').toLowerCase()
            assert normalised_json.containsKey(normalised_name)
            def db_version = normalised_json[normalised_name]
            def db_dir_parts = app_dirs[normalised_name]
            def db_dir = db_dir_parts[0]
            def db_subdir = db_dir_parts[1]
            Path path = Paths.get("${datadir}/${db_dir}/${db_version}/${db_subdir}")
            if (Files.exists(path)) {
                with_data.add( [ normalised_name, db_version, path.toString() ])
            } else {
                without_data.add( [ normalised_name, db_dir, db_version, to_download.contains(db_dir), path.toString() ] )
                to_download.add( db_dir )
            }
        }
    }
}

process VALIDATE_DATA {
    label    'tiny'
    executor 'local'
    cache false  // Stops the esotericsoftware.kryo.serializers warning

    input:
    val list_databases

    output:
    val map_databases

    exec:
    map_databases = list_databases.collectEntries { name, version, dirpath ->
        [(name): [version: version, dirpath: dirpath]]
    } 
}
