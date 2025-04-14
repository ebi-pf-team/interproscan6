import java.io.File
import java.nio.file.*
import groovy.json.JsonSlurper

process DOWNLOAD {
    maxForks 1
    label    'local', 'ips6_container'

    input:
    tuple val(name), val(version)
    val iprscan_version
    val outdir

    output:
    val true

    script:
    """
    cd ${outdir}
    curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${name}/${name}-${version}.tar.gz
    curl -OJ ${InterProScan.FTP_URL}/${iprscan_version}/${name}/${name}-${version}.tar.gz.md5
    md5sum -c ${name}-${version}.tar.gz.md5 || { echo "Error: MD5 checksum failed" >&2; exit 1; }
    tar -zxf ${name}-${version}.tar.gz
    rm ${name}-${version}.tar.gz*
    chmod 777 -R ${name}
    """
}

process VALIDATE_APP_DATA {
    input:
    val ready
    val json_file
    val databases
    val app_dirs
    val datadir

    output:
    path "missing_databases.json"

    exec:
    File file = new File(json_file)
    def json = new JsonSlurper().parse(file)
    def normalised_json = [:]
    json.each { key, value ->
        normalised_json[key.replaceAll(/[\s\-]+/, '').toLowerCase()] = value
    }

    databases_with_version = [] as Set
    databases.each { db_name ->
        if (app_dirs.containsKey(db_name)) {
            def normalised_name = db_name.replaceAll(/[\s\-]+/, '').toLowerCase()
            assert normalised_json.containsKey(normalised_name)
            def db_version = normalised_json[normalised_name]
            def db_dir_parts = app_dirs[normalised_name]
            def db_dir = db_dir_parts[0]
            def db_subdir = db_dir_parts[1]
            Path path = Paths.get("${datadir}/${db_dir}/${db_version}/${db_subdir}")
            if (!Files.exists(path)) {
                databases_with_version.add( [db_dir, db_version] )
            }
        }
    }

    // Pass the output to a file to avoid the esotericsoftware.kryo.serializers warning
    def outputFile = task.workDir.resolve("missing_databases.json")
    outputFile.text = groovy.json.JsonOutput.toJson(databases_with_version)
    return outputFile.toString()
}

process GET_DB_RELEASES {
    // Run as a process to ensure it runs after all downloading is complete
    input:
    val db_json_file
    val xref_dir
    val xref_config
    val add_goterms
    val add_pathways
    val ready

    output:
    path "databases_with_version.json"

    exec:
    // Before getting all the dbs (inc. interpro) version numbers check we have all the interpro data
    // member data was checked in PREPARE_DOWNLOADS
    error = InterProScan.validateXrefFiles(
        xref_dir,
        xref_config,
        add_goterms,
        add_pathways
    )
    if (error) {
        log.error error
        exit 1
    }

    JsonSlurper jsonSlurper = new JsonSlurper()
    def databaseJson = new File(db_json_file.toString())
    databases_with_version = jsonSlurper.parse(databaseJson)
    databases_with_version = databases_with_version.collectEntries { appName, versionNum ->
        [(appName.toLowerCase()): versionNum]
    }

    // Pass the output to a file to avoid the esotericsoftware.kryo.serializers warning
    def outputFile = task.workDir.resolve("databases_with_version.json")
    outputFile.text = groovy.json.JsonOutput.toJson(databases_with_version)
    return outputFile.toString()
}