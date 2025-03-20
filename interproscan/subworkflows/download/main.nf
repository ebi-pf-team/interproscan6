workflow DOWNLOAD_INTERPRO {
    take:
    applications           // list of applications to run
    appsConfig             // map of applications
    datadir                // path to the interpro data directory

    main:
    /* Get the iprscan version as this determines which InterPro versions are compatiable
    Only the major and minor versions are reported on the ftp. Changing the minor version
    indicates incompatibility */
    def iprscan_version = workflow.manifest.version.split("\\.")[0..1].join(".")
    def ftpUrl = InterPro.sanitizeURL("${params.ftpServer.url}/${iprscan_version}")

    // Check if /xrefs/databases.json exists
    def databaseJson = "$datadir/${params.xRefsConfig.entries}"
    if (!InterProScan.resolveFile(databaseJson)) {
        // No <datadir>/xrefs/databases.json so download a fresh InterPro dataset
        buildDataDir(datadir)
        // [1] Get the latest, compatible InterPro version for this Iprscan version
        def response = InterPro.httpRequest(ftpUrl, null, 0, true, log, false)
        def latestRelease = extractHrefs(response).max()
        // [2] Process each application in turn
    } else {
        // Download data for the existing InterPro dataset
        // [1] Get current InterPro release version
        def interpro_version = InterPro.getDatabaseVersion("InterPro", datadir).toFloat()

        // [2] Get list of all compatiable InterPro releases
        def response = InterPro.httpRequest(ftpUrl, null, 0, true, log, false)
        def releases = extractHrefs(response)

        // [3] Check if releases are compatiable
        if (!releases.contains(interpro_version)) {
            log.error "InterProScan version ${workflow.manifest.version} is not compatiable with InterPro release ${interpro_version}"
            exit 1
        } else {
            println "nothing to do here yet"
            // [4] Check per member db which data needs to be downloaded
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