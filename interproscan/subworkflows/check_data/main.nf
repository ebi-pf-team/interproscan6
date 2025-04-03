workflow CHECK_DATA {
    /* "Ensure the required data is present in the data directory and member database files.
    If any are missing, terminate the process. */
    take:
    apps              // member databases to run
    interproRelease   // str, interproRelease if --download enabled else an empty string

    main:
    def _datadir = ""
    def _interproDir = ""
    def _interproRelease = ""
    def _iprScanVersion = ""
    def _memberDbReleases = [:]

    // Validate the data dir and application data files if needed by any members
    // e.g. mobidblite and coils do no need additional data files
    need_data = apps.any { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }

    if (need_data) {
        // Check if there is a data directory
        // if --datadir is called and no path provided it can revert to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            log.error error
            exit 1
        }

        // Get the number of the InterPro release to use
        if (!params.offline) {  // Online - check the interpro-iprscan compatibility
            if (!interproRelease) {
               (_interproRelease, error) = InterPro.getInterproRelease(_datadir, params.interpro.toString())
                if (error) {
                    log.error error
                    exit 1
                }
            } else {
                _interproRelease = interproRelease
            }

            (noConnWarning, latestReleaseWarning, compatError) = InterPro.checkCompatibility(_iprScanVersion, _interproRelease, log)
            if (compatError) { // Not compatible, terminate
                log.error compatError
                exit 1
            } else if (noConnWarning) {  // Could not connect to the ftp, so continuing assuming it's compatible
                log.warn noConnWarning
            } else if (latestReleaseWarning) {  // Let the user know a later InterPro release is available
                log.warn latestReleaseWarning
            }
        } else {  // Offline - do not check for interpro-iprscan compatibility
            (_interproRelease, error) = InterPro.getInterproRelease(_datadir, params.interpro.toString())
            if (error) {
                log.error error
                exit 1
            }
        }
        println "# InterPro $_interproRelease"  // adds to the pipeline intro lines in main.nf

        // Validate the interpro and Xref data files
        // Validate these first as we need the interpro files for validateAppData
        error = InterProScan.validateXrefFiles(_datadir, _interproRelease, params.xRefsConfig, params.goterms, params.pathways)
        if (error) {
            log.error error
            exit 1
        }

        /* Load the database.json file and set all keys to lowercase to match applications.config
        Don't worry about checking it exists, this was done in InterProScan.validateXrefFiles() */
        _memberDbReleases = InterPro.getMemberDbReleases(params.xRefsConfig, _interproRelease, _datadir)

        // Validate the selected member databases data files and dirs
        error = InterProScan.validateAppData(apps, _datadir, params.appsConfig, _memberDbReleases)
        if (error) {
            log.error error
            exit 1
        }
    } // end of (need_data)

    datadir = _datadir
    finalInterproRelease = _interproRelease
    memberDbReleases = _memberDbReleases

    emit:
    datadir                // str: path to data directory
    finalInterproRelease   // str: interpro db version number
    memberDbReleases       // map: [db: version number]
}