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

    // Download data - if selected

    // Validate the data dir and application data files if needed by any members
    // e.g. mobidblite and coils do no need additional data files
    need_data = apps.any { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }
    def _datadir = ""
    if (need_data) {
        // If --datadir is called and no path is given it converts to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            log.error error
            exit 1
        }
        error = InterProScan.validateAppData(apps, _datadir, params.appsConfig)
        if (error) {
            log.error error
            exit 1
        }
        error = InterProScan.validateXrefFiles(_datadir, params.xRefsConfig, params.goterms, params.pathways)
        if (error) {
            log.error error
            exit 1
        }
    }

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

    // Sequences validation
    (error, numSequences) = FastaFile.validate(params.input, params.nucleic, params.appsConfig, apps)
    if (error) {
        log.error error
        exit 1
    } else if (!numSequences) {
        log.error "No FASTA sequences found in ${params.input}"
        exit 1
    }

    def _matchesApiUrl = null
    if (params.offline && params.matchesApiUrl != null) {
        log.error "--offline and --matches-api-url are mutually exlusive"
        exit 1
    } else if (params.offline) {
        _matchesApiUrl = null
    } else {
        _matchesApiUrl = InterPro.getMatchesApiUrl(
            params.matchesApiUrl, params.lookupService.url, _datadir, workflow.manifest, log
        )
    }

    matchesApiUrl = _matchesApiUrl
    datadir = _datadir

    emit:
    fasta         // str: path to input fasta file
    datadir       // str: path to data directory
    apps          // list: list of application to
    outdir        // str: path to output directory
    formats       // set<String>: output file formats
    signalpMode   // str: Models to be used with SignalP
    matchesApiUrl // str|null: URL of Matches API to query
}