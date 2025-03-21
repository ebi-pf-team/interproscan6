import groovy.json.JsonSlurper

include { DOWNLOAD_INTERPRO } from "../download"
include { getMatchesApiUrl  } from "../../modules/lookup"

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
    def download = false
    if (need_data) {
        // [1] Check the data dir exists
        // If --datadir is called and no path is given it converts to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (_datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!_datadir) {
            if (params.download && !params.datadir) {
                log.error "--datadir <DATA-DIR> is mandatory when using --download"
                exit 1
            } else if (params.download) {
                download = true
            } else { // datadir is missing, it's needed and --download is not enabled
                log.error error
                exit 1
            }
        }

        // [2] Check datadir/versions.json exists or download the file if able
        // This file is essential is it ensures the InterPro and Iprscan versions are compatible
        def _versionsFile = _datadir.resolve("versions.json").toString())
        if (!InterProScan.resolveFile(_versionsFile)) {
            if (params.download) {
                // download versions.json /* TODO: Write download version.json function */
                InterPro.downloadVersions(_versionsFile)
            } else {
                // datadir/versions.json is missing, it's needed and --download is not enabled
                log.error ("Could not find 'versions.json' in the data directory ${_datadir}.\n" +
                           "Please ensure that the data dir is correctly populated or use --download")
                exit 1
            }
        }

        // [3] Get the iprscan version, which determines the compatible interpro releases
        // Only the major and minor versions are needed
        def _iprScanVersion = workflow.manifest.version.split("\\.")[0..1].join(".")

        // [4] Get the version of InterPro to be used
        def _interproVersion = null
        def _interproDir = ""
        def _versionsMaps = new JsonSlurper().parse(new File(_versionsFile.toString())) // E.g. [6.0: [102.0, 103.0, 104.0]]
        if (!_versionsMaps.containsKey(_iprScanVersion)) {
            log.error "InterProScan version ${_iprScanVersion} is not listed in ${_versions}"
            exit 1
        } else {
            if (params.interpro == "latest") {
                // Get the list of compatible InterPro releases
                def _compatibleVersions = _versionsMaps[_iprScanVersion]*.toFloat()
                // Check which InterPro releases are already in the datadir, if any
                def _interproDir = new File(_datadir.resolve("interpro").toString())
                if (_interproDir.exists() && _interproDir.isDirectory()) {
                    def _dirs = _interproDir.listFiles()
                                .findAll { it.isDirectory() }
                                .collect { it.name.toFloat }
                                .findAll { it in _compatibleVersions }
                    if (!_dirs) {
                        // No compatible InterPro releases are stored in the data dir
                        if (params.download) {
                            download = true
                        } else {
                            // No compatible InterPro releases and --download is not enabled
                            log.error ("No compatible InterPro releases are stored in the data dir ${_datadir}.\n" +
                                       "Compatible releases: ${_compatibleVersions.join(', ')}\n" +
                                       "Please ensure that the data dir is correctly populated or use --download")
                            exit 1
                        }
                    } else {
                        _latestDir = _dirs.max()
                        if (_latestDir < _compatibleVersions.max()) {
                            if (params.download) {
                                _interproVersion = _compatibleVersions.max()
                                download = true
                            } else {
                                _interproVersion = _dirs.max()
                                log.info("A newer compatible InterPro release (${_compatibleVersions.max()}) is available.\n" +
                                         "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release.")
                            }
                        } else {
                            _interproVersion = _compatibleVersions.max()
                        }
                    }
                } else {
                    // _interproDir does not exist
                    if (params.download) {
                        _interproVersion = _compatibleVersions
                        download = true
                    } else {
                        log.error ("No 'interpro' directory was found in the data directory ${_datadir}\n" +
                                   "Please ensure that the data dir is correctly populated or use --download")
                        )
                        exit 1
                    }
                }
            } else {
                // Using the user selected InterPro release. But first check it is compatible with the Iprscan version
                if (_versionsMaps[_iprScanVersion].contains(params.interpro.toString())) {
                    _interproVersion = params.interpro.toString()
                } else {
                    log.error ("Interpro version ${params.interpro} is not compatible with InterProScan " +
                               "version ${_iprScanVersion}")
                    exit 1
                }
            }
        }
        println "# InterPro: $_interproVersion"

        // [5] Download data if needed and enabled
        if (download) {
            DOWNLOAD_INTERPRO(apps, params.appsConfig, params.datadir, _interproVersion, _iprScanVersion)
        }

        // [6] Validate the selected member databases data files and dirs
        error = InterProScan.validateAppData(apps, _datadir, params.appsConfig)
        if (error) {
            log.error error
            exit 1
        }

        // [7] Validate the interpro Xref data files
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
        log.error "--offline and --matches-api-url are mutually exclusive"
        exit 1
    } else if (params.offline) {
        _matchesApiUrl = null
    } else {
        _matchesApiUrl = getMatchesApiUrl(
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