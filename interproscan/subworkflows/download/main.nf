include { DOWNLOAD as DOWNLOAD_INTERPRO
          DOWNLOAD as DOWNLOAD_MEMBERDB
          CHECK_APP_DATA                } from "../../modules/download"
// Load with different names to make the terminal output clearer for the user

import java.io.File

workflow DOWNLOAD_DATA {
    // This subworkflow runs when the --download flag is used.

    take:
    iprScanVersion         // str, major and minor iprscan release number
    applications           // list of applications to run
    appsConfig             // map of applications

    main:
    def baseURL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/${iprScanVersion}"
    def _datadir = ""
    def dirPath = ""
    def downloadInterPro = false
    def downloadURL = ""
    def interproReleasesMap = [:] // Downloaded from the ftp = [interpro: [release(s)]]
    def savePath = ""
    def toDownload = [] // names of member dbs to download
    def _interproRelease = null

    // [1] Identify the apps that need data to run
    _apps = applications.findAll { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }
    println "DEBUG: [1]"

    // [2] Get the list of InterPro releases that are compatible with this InterProScan release
    def _releasesURL = "${baseURL}/versions.json"
    interproReleasesMap = InterPro.httpRequest(_releasesURL, null, 4, true, log)
    if (!interproReleasesMap) {
        log.error "Failed to retrieve the InterPro release 'versions.json' file from the FTP"
        exit 1
    }
    def compatibleReleases = interproReleasesMap["interpro"]*.toFloat()
    println "DEBUG: [2]"

    // [3] Get the InterPro release to be used
    if (params.interpro == "latest") {
        _interproRelease = compatibleReleases.max()
        println "DEBUG: [3a]"
    } else {
        // Using the user selected InterPro release, but first check it is compatible with this iprScan release
         _interproRelease = params.interpro.toString()
         if (!compatibleReleases.contains(_interproRelease)) {
             log.error "The InterPro release ${_interproRelease} is not compatible with InterProScan version ${iprScanVersion}"
             exit 1
         } else if (params.interpro.toFloat() < compatibleReleases.max()) {
            log.info "A later compatible InterPro release '${compatibleReleases[iprscan].max()}' is availble.\n"+
                     "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release"
         }
         println "DEBUG: [3b]"
    }
    println "DEBUG: [3c]"

    // [4] Do we need to download InterPro data?
    // [4a] Check if the datadir exists
    // If --datadir is called and no path is given, it can convert to a boolean
    dirPath = params.datadir instanceof Boolean ? null : params.datadir
    (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
    if (!_datadir) {
        // Will need to download everything
        downloadInterPro = true
        toDownload = _apps  // download data for all apps that need data
        _datadir = new File(dirPath)
        _datadir.mkdirs()
        println "DEBUG: [4a]"
    } else {
        // [4b] Check if we need to specifically download the InterPro data
        def interproDir = new File(_datadir.resolve("interpro/$_interproRelease").toString())
        if (!interproDir.exists()) {
            downloadInterPro = true
        } else {
            // check all the files are present in <datadir>/interpro/<version>/*
            error = InterProScan.validateXrefFiles(_datadir, _interproRelease.toString(), params.xRefsConfig, params.goterms, params.pathways)
            if (error) {
                downloadInterPro = true
            }
        }
        println "DEBUG: [4b]"
    }
    println "DEBUG: [4c]"

    // [5] Download InterPro
    println "DEBUG: DOWNLOAD: $downloadInterPro"
    if (downloadInterPro) {
        log.info "Downloading InterPro release $_interproRelease"
        DOWNLOAD_INTERPRO(["interpro", _interproRelease, baseURL, _datadir], _interproRelease).collect()
        // [6] Compile params for downloading data for the member dbs with missing data
        CHECK_APP_DATA(DOWNLOAD_INTERPRO.out, baseURL, _interproRelease, _apps, params.appsConfig, params.xRefsConfig)
    } else {
        // [6] Compile params for downloading data for the member dbs with missing data
        CHECK_APP_DATA(_datadir, baseURL, _interproRelease, _apps, params.appsConfig, params.xRefsConfig)
        // Do I need to process download_params?
    }
    println "DEBUG: [5] & [6]"
    
    download_params = CHECK_APP_DATA.out.flatMap()

    // [7] Download member database data if needed
    if (download_params) {
        interproRelease = DOWNLOAD_MEMBERDB(download_params, _interproRelease)
    } else {
        interproRelease = _interproRelease
    }
    println "DEBUG: [7]"
    log.info "DEBUG: interproRelease - $interproRelease"
    log.error "DEBUG: Got to end of DOWNLOAD_DATA"
//     exit 2

    emit:
    interproRelease       // str: InterPro release version number
}
