workflow DOWNLOAD_INTERPRO {
    // This subworkflow runs when the --download flag is used.

    take:
    applications           // list of applications to run
    appsConfig             // map of applications

    main:
    def baseUrl = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"
    def _datadir = ""
    def dirPath = ""
    def dbReleasesMap = [:] // Downloaded from the ftp = [db: [release(s)]]
    def downloadURL = ""
    def dbReleasesMap = [:]
    def savePath = ""
    def toDownload = [] // names of dbs to download
    def _interproRelease = null

    // Apps that have data to download
    _apps = applications.findAll { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }

    // Establish which data needs to be downloaded
    while (toDownload.isEmpty()) {
        // [1] Check the data dir exists. Do we need to build one?
        // If --datadir is called and no path is given it converts to a boolean
        dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            buildDataDir(_datadir)
            toDownload << "interpro"
            downloadApps = _apps  // download data for all apps that need data
            break  // go straight to downloading everything
        }

        // [2] Get the iprscan version (only major and minor version needed), to determine the compatible interpro releases
        def _iprScanVersion = workflow.manifest.version.split("\\.")[0..1].join(".")

        // [3] Get the list of InterPro releases that are compatible with this InterProScan release
        def _releasesURL = "${baseUrl}/${_iprScanVersion}/versions.json"
        def dbReleasesMap = InterPro.httpRequest(_releasesURL, null, 0, true, log)
        def compatibleReleases = dbReleasesMap["interpro"]*.toFloat()

        // [4] Get the InterPro release to be used
        if (params.interpro == "latest") {
            _interproRelease = compatibleReleases.max()
        } else {
            // Using the user selected InterPro release, but first check it is compatible with this iprScan release
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

        // [5] Check if we need to download the InterPro and XREF files
        def interproDir = new File(_datadir.resolve("interpro/$_interproRelease").toString())
        if (!interproDir.exists()) {
            toDownload << "interpro"
        } else {
            error = InterProScan.validateXrefFiles(_datadir, params.xRefsConfig, params.goterms, params.pathways)
            if (error) {
                toDownload << "interpro"
            }
        }

        // [6] Check if we need to download any member db data. validateAppData returns a set of apps with missing data
        downloadApps.addAll(InterProScan.validateAppData(_apps, _datadir, params.appsConfig, returnSet=true))
        break
    }

    def dbRelease = ""
    toDownload.each { dbName ->
        dbRelease = (dbName == "interpro") ? _interproRelease : dbReleasesMap[dbName]
        downloadURL = "${baseUrl}/$dbName/$adbName-$dbRelease.tar.gz"
        savePath = _datadir.resolve("$adbName/$adbName-$dbRelease.tar.gz")
        dirPath = _data.resolve("$adbName/$dbRelease")
        error = InterPro.downloadFile(downloadURL, savePath)
        if (error) {
            log.error error
            exit 1
        }
        error = InterPro.extractTarFile(savePath, dirPath)
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
