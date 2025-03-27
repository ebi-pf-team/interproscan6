workflow DOWNLOAD_INTERPRO {
    // This subworkflow runs when the --download flag is used.

    take:
    applications           // list of applications to run
    appsConfig             // map of applications

    main:
    def _datadir = ""
    def _downloadInterpro = false
    def _downloadApps = false
    def _interproReleases = [:]
    def _interproRelease = null
    def baseUrl = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"

    // Apps that have data to download
    _apps = applications.findAll { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }

    // Establish which data needs to be downloaded
    while (!_downloadInterpro && !_downloadApps) {
        // [1] Check the data dir exists
        // If --datadir is called and no path is given it converts to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            if (!params.datadir) {
                log.error "--datadir <DATA-DIR> is mandatory when using --download"
                exit 1
            } else { // datadir is missing, it's needed and --download is enabled
                _downloadInterpro = true
                _downloadApps = true
                break  // go straight to downloading everything
            }
        }

        // [2] Get the iprscan version, which determines the compatible interpro releases
        // Only the major and minor versions are needed
        def _iprScanVersion = workflow.manifest.version.split("\\.")[0..1].join(".")

        // [3] Get the list of InterPro releases that are compatible with this InterProScan release
        def _releasesURL = "${baseUrl}/${_iprScanVersion}/versions.json"
        def _interproReleases = InterPro.httpRequest(_releasesURL, null, 0, true, log)
        def _compatibleReleases = _interproReleases[_iprScanVersion]*.toFloat()

        // [4] Get the InterPro release to be used
        if (params.interpro == "latest") {
            _interproRelease = _compatibleReleases.max()
        } else {
            // Using the user selected InterPro release. But first check it is compatible with this iprScan release
             _interproRelease = params.interpro.toString()
             if (!_compatibleReleases.contains(_interproRelease)) {
                 log.error "The InterPro release ${_interproRelease} is not compatible with InterProScan version ${_iprScanVersion}"
                 exit 1
             }
        }
        println "# InterPro: $__interproRelease" // Adds to the summary statements in the main workflow

        // [5] Check if the InterPro (entries and XREF) data is available locally
        def _interproDir = new File(_datadir.resolve("interpro").resolve(_interproRelease).toString())
        if (!_interproDir.exists() || !_interproDir.isDirectory()) {
            _downloadInterpro = true
        }

        // [6] Check if we need to download any member database data
        // InterProScan.validateAppData() only returns something if there is missing data
        if (InterProScan.validateAppData(_apps, _datadir, params.appsConfig, true)) {
            _downloadApps = true
        }

        break
    }

    if (_downloadInterpro) {
        // Download the InterPro (entries+xref) data
    }

    if (_downloadApps) {
        // Download member database data
    }

    interproRelease = _interproRelease

    emit:
    interproRelease       // str: InterPro release version number
}

def buildDataDir(String datadir) {
    // Build a directory to store data
    return
}

def validateAppData(List<String> apps, String baseUrl) {
    apps.each { String app ->
        // Compile path to app datadir
        def appDataDir = "$datadir/$app"

        def downloadURL = ""
        def cmd = []
        def process = cmd.execute()
        process.waitFor()

        // Check the md5 hash to ensure file was correctly download

    }
}

// Get paths for files that need downloading
// If files are present check the md5 hash
// If files are not present or md5 hash does not match