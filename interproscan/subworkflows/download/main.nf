import java.net.URL
import java.net.MalformedURLException
import jave.io.File
import java.io.FileOutputStream
import java.io.IOException
import java.io.InputStream
import java.io.OutputStream

workflow DOWNLOAD_INTERPRO {
    // This subworkflow runs when the --download flag is used.

    take:
    applications           // list of applications to run
    appsConfig             // map of applications

    main:
    def _datadir = ""
    def dbReleasesMap = [:] // Downloaded from the ftp = [db: [release(s)]]
    def downloadInterpro = false
    def downloadApps = []
    def dbReleasesMap = [:]
    def _interproRelease = null
    def baseUrl = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"
    def downloadURL = ""
    def savePath = ""

    // Apps that have data to download
    _apps = applications.findAll { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }

    // Establish which data needs to be downloaded
    while (!downloadInterpro && !downloadApps) {
        // [1] Check the data dir exists. Do we need to build one?
        // If --datadir is called and no path is given it converts to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            buildDataDir(_datadir)
            downloadInterpro = true
            downloadApps = _apps  // download data for all apps that need data
            break  // go straight to downloading everything
        }

        // [2] Get the iprscan version, which determines the compatible interpro releases
        // Only the major and minor versions are needed
        def _iprScanVersion = workflow.manifest.version.split("\\.")[0..1].join(".")

        // [3] Get the list of InterPro releases that are compatible with this InterProScan release
        def _releasesURL = "${baseUrl}/${_iprScanVersion}/versions.json"
        def dbReleasesMap = InterPro.httpRequest(_releasesURL, null, 0, true, log)
        def compatibleReleases = dbReleasesMap["interpro"]*.toFloat()

        // [4] Get the InterPro release to be used
        if (params.interpro == "latest") {
            _interproRelease = compatibleReleases.max()
        } else {
            // Using the user selected InterPro release. But first check it is compatible with this iprScan release
             _interproRelease = params.interpro.toString()
             if (!compatibleReleases.contains(_interproRelease)) {
                 log.error "The InterPro release ${_interproRelease} is not compatible with InterProScan version ${_iprScanVersion}"
                 exit 1
             } else if (params.interpro.toFloat() < compatibleReleases.max()) {
                log.info "A later compatible InterPro release '${compatibleReleases[iprscan].max()}' is availble.\n"+
                         "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release"
             }
        }
        println "# InterPro: $_interproRelease" // Adds to the summary statements in the main workflow

        // [7] Check if we need to download the correct InterPro and XREF files
        // TODO: check the MD5 check
        def xrefMd5s = InterProScan.validateXrefFiles(datadir, params.xRefsConfig, params.goterms, params.pathways, returnList=true)
        xrefMd5s.each {String xrefFile ->
            if (!xrefMd5s[xrefFile]) {
                // File not found
                downloadInterpro = true
            } else {
                // Check the md5 is correct
                 (ftpMd5, error) = InterPro.getFtpMd5("interpro", _interproRelease.toString())
                 if (error) {
                    log.error error
                    exit 1
                 } else if (xrefMd5s[xrefFile] != ftpMd5) {
                    downloadInterpro = true
                 }
            }
        }

        // [6] Check if we need to download any member database data
        // InterProScan.validateAppData() only returns something if there is missing data
        // TODO: Add a check md5 sum check
        downloadApps = InterProScan.validateAppData(_apps, _datadir, params.appsConfig, returnSet=true)
        break
    }

    if (downloadInterpro) {
        downloadURL = "${baseUrl}/interpro/interpro-${_interproRelease}.tar.gz"
        savePath = _datadir.resolve("interpro/interpro-${_interproRelease}.tar.gz")
        error = InterPro.downloadFile(downloadURL, savePath)
        if (error) {
            log.error error
            exit 1
        }
    }

    def dbRelease = ""
    downloadApps.each { app ->
        dbRelease = dbReleasesMap[app]
        downloadURL = "${baseUrl}/app/$app-$dbRelease.tar.gz"
        savePath = _datadir.resolve("$app/$app-${dbRelease}.tar.gz")
        error = InterPro.downloadFile(downloadURL, savePath)
        if (error) {
            log.error error
            exit 1
        }
    }

    interproRelease = _interproRelease

    emit:
    interproRelease       // str: InterPro release version number
}

def buildDataDir(String dirPath) {
    def dir = new File(dirPath)
    dir.mkdirs()
}
