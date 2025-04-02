import java.io.File
import java.nio.file.*

process DOWNLOAD {
    maxForks 1
    label    'local', 'ips6_container'

    input:
    tuple val(database), val(dbVersion), val(baseURL), val(datadir)
    // DATADIR must be a path object so that we can build the absolute path

    output:
    val datadir

    script:
    def d = datadir
    def v = dbVersion
    def db = database
    println "DEBUG: ${d}: ${d.getClass()}. $db. $v"
    def tarURL = "${baseURL}/${database}/${database}-${dbVersion}.tar.gz"
    def tarSavePath = "${database}-${dbVersion}.tar.gz"
    def md5URL = "${tarURL}.md5"
    def md5SavePath = "${tarSavePath}.md5"
    // Convert the final extract path to an absolute path to ensure files are written to the datadir
    def downloadDirPath = datadir.getAbsolutePath()
    println "DEBUG: ==\n${tarURL}\n$md5URL\n$tarSavePath\n$md5SavePath\n$downloadDirPath\n=="
    """
    mkdir -p ${downloadDirPath}  # ensure that the datadir exists
    wget -O ${tarSavePath} "${tarURL}"
    wget -O ${md5SavePath} "${md5URL}"

    DOWNLOADED_MD5=\$(md5sum ${tarSavePath} | awk '{ print \$1 }')
    PRECOMPUTED_MD5=\$(cat ${md5SavePath} | awk '{ print \$1 }')

    if [ "\$DOWNLOADED_MD5" != "\$PRECOMPUTED_MD5" ]; then
        echo "Error: MD5 checksum mismatch for ${tarURL}"
        exit 1
    else
        # Extract to produce datadir/member/version/*
        tar -xvf ${tarSavePath} -C ${downloadDirPath}
        rm ${tarSavePath}
    fi
    """
}

process CHECK_APP_DATA {
    // This process is used to force the sequential operation between download InterPro and member db data
    maxForks 1
    label    'local'

    input:
    val datadir          // absolute path to --datadir in an ArrayBag
    val baseURL
    val interproRelease
    val applications
    val appsConfig
    val xRefsConfig

    output:
    val downloadParams

    exec:
    def d = datadir
    println "datadir -- ${d.getClass()} -- $d -- ${d[0]}"
    def datadirPath = new File(datadir[0].toString())
    println "datadirPath --${datadirPath.getClass()} -- $datadirPath"
    def memberDbReleases = InterPro.getMemberDbReleases(xRefsConfig, interproRelease.toString(), datadirPath)
    def appsToDownload = InterProScan.validateAppData(applications as List, datadir[0] as Path, appsConfig, memberDbReleases, returnSet=true)

    downloadParams = appsToDownload.collect { app ->
        def appVersion = memberDbReleases[app]
        return [app, appVersion, baseURL, datadir[0]]
    }
    println "DEBUG: $downloadParams"
    return downloadParams
}