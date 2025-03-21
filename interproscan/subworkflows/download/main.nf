workflow DOWNLOAD_INTERPRO {
    take:
    applications           // list of applications to run
    appsConfig             // map of applications
    datadir                // path to the interpro data directory

    main:
    /* Get the iprscan version as this determines which InterPro versions are compatible
    Only the major and minor versions are reported on the ftp */
    def iprscan_version = workflow.manifest.version.split("\\.")[0..1].join(".")
    def ftpUrl = InterPro.sanitizeURL("${params.ftpServer.url}/${iprscan_version}")

    // Apps that have data to download
    _apps = applications.findAll { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }

    // Check if /xrefs/databases.json exists
    def databaseJson = "$datadir/${params.xRefsConfig.entries}"
    if (!InterProScan.resolveFile(databaseJson)) {
        // No <datadir>/xrefs/databases.json so download a fresh InterPro dataset
        // [1] Build the datadir (if needed)
        buildDataDir(datadir)

        // [2] Get the latest, compatible InterPro version for this Iprscan version
        def response = InterPro.httpRequest(ftpUrl, null, 0, true, log, false)
        def interpro_version = extractHrefs(response).max()
        log.info "Downloading the latest compatible InterPro release ($interpro_version) to $datadir"

        // [3] Compile the parent download URL
        def parentURL = "$ftpURL/$interpro_version"

        // [4] Process each application in turn

    } else {
        // Download data for the existing InterPro dataset

        // [1] Get current InterPro release version
        def interpro_version = InterPro.getDatabaseVersion("InterPro", datadir).toFloat()

        // [2] Get list of all compatiable InterPro releases
        def response = InterPro.httpRequest(ftpUrl, null, 0, true, log, false)
        def releases = extractHrefs(response)

        // [3] Check if releases are compatiable
        if (!releases.contains(interpro_version)) {
            log.error "InterProScan version ${workflow.manifest.version} is not compatible with InterPro release ${interpro_version}"
            exit 1
        } else {
            log.info "Downloading InterPro release ${interpro_version} to ${datadir}"
            // [4] Compile the parent download URL
            def parentURL = "$ftpURL/$interpro_version"

            // [5] Check per member db which data needs to be downloaded

        }
    }
}

def buildDataDir(String datadir) {
    // Build a directory to store data
    return
}

List<String> extractHrefs(String htmlContent) {
    // Extract the hrefs from a ftp page, e.g. from https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0/
    def hrefs = []
    def matcher = htmlContent =~ /<img [^>]*alt="\[DIR\]"[^>]*>.*?<a href="([^"]+)\/"/
    matcher.each { match ->
        try {
            hrefs << match[1].toFloat()
        } catch (NumberFormatException e) {
            // Skip if conversion fails. In case other dirs are added in the future
        }
    }
    return hrefs
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