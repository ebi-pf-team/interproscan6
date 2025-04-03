include { DOWNLOAD_DATA } from "../download"
include { getMatchesApiUrl } from "../../modules/lookup"

workflow INIT_PIPELINE {
    main:
    // Params validation
    InterProScan.validateParams(params, log)

    // Check the input
    fasta = InterProScan.resolveFile(params.input)
    if (!fasta) {
        log.error "No such file: ${params.input}"
        exit 1
    }

    // Applications validation
    (apps, error) = InterProScan.validateApplications(params.applications, params.appsConfig)
    if (!apps) {
        log.error error
        exit 1
    }

    if (params.download && params.offline) {
        log.error "--download and --offline are mutually exclusive"
        exit 1
    }

    // Validate the data dir and application data files if needed by any members
    // e.g. mobidblite and coils do no need additional data files
    need_data = apps.any { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }
    def _datadir = ""
    def _interproDir = ""
    def _interproRelease = ""
    def _iprScanVersion = ""
    def _memberDbReleases = [:]
    if (need_data) {
        // Get the iprscan version (only major and minor version needed), to determine the compatible interpro releases
        _iprScanVersion = workflow.manifest.version.split("\\.")[0..1].join(".")

        if (params.download) {
            if (!params.datadir) {
                log.error "--datadir <DATA-DIR> is mandatory when using --download"
                exit 1
            }
            println "*** Downloading "
            DOWNLOAD_DATA(_iprScanVersion, apps, params.appsConfig)
            println "******** Downloading"
            _interproRelease = DOWNLOAD_DATA.out.interproRelease.val
        } else {
            // Check if there is a data directory
            // If --datadir is called and no path is given it converts to a boolean
            def dirPath = params.datadir instanceof Boolean ? null : params.datadir
            (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
            if (!_datadir) {
                log.error error
                exit 1
            }

            // Get the number of the InterPro release to use
            if (!params.offline) {  // Online - check the interpro-iprscan compatibility
                (_interproRelease, error) = InterPro.getInterproRelease(_datadir, params.interpro.toString())
                if (error) {
                    log.error error
                    exit 1
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

    // Check valid output file formats were provided
    (formats, error) = InterProScan.validateFormats(params.formats)
    if (error) {
        log.error error
        exit 1
    }

    // Build output dir if needed
    (outdir, error) = InterProScan.resolveDirectory(params.outdir, false, true)
    if (!outdir) {
        log.error error
        exit 1
    }

    // SignalP mode validation
    (signalpMode, error) = InterProScan.validateSignalpMode(params.signalpMode)
    if (!signalpMode) {
        log.error error
        exit 1
    }

    def _matchesApiUrl = null
    if (params.offline && params.matchesApiUrl != null) {
        log.error "--offline and --matches-api-url are mutually exclusive"
        exit 1
    } else if (params.offline) {
        _matchesApiUrl = null
    } else {
        _matchesApiUrl = getMatchesApiUrl(
            params.matchesApiUrl, params.lookupService.url, _interproRelease, workflow.manifest, log
        )
    }

    matchesApiUrl = _matchesApiUrl
    datadir = _datadir
    interproRelease = _interproRelease
    memberDbReleases = _memberDbReleases

    emit:
    fasta            // str: path to input fasta file
    datadir          // str: path to data directory
    interproRelease  // str: interpro db version number
    memberDbReleases // map: [db: version number]
    apps             // list: list of application to
    outdir           // str: path to output directory
    formats          // set<String>: output file formats
    signalpMode      // str: Models to be used with SignalP
    matchesApiUrl    // str|null: URL of Matches API to query
}
