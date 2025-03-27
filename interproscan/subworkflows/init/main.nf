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
    def _interproDir = ""
    def _interproRelease = ""
    if (need_data) {
        if (params.download) {
            DOWNLOAD_INTERPRO(apps, params.appsConfig)
            _interproRelease = DOWNLOAD_INTERPRO.out.interproRelease.val
        } else if (!params.offline) {
            (_interproDir, error) = InterProScan.validateInterproDir(_datadir)
            if (error) {  // no datadir/interpro dir found
                log.error error
                exit 1
            }
            _interproRelease = InterPro.getInterproRelease(_datadir, params.interpro.toString())

            // Check compatibility

        } else {
            // Working offline so do not check for interpro-iprscan compatibility
            (_interproDir, error) = InterProScan.validateInterproDir(_datadir)
            if (error) {  // no datadir/interpro dir found
                log.error error
                exit 1
            }
            _interproRelease = InterPro.getInterproRelease(_datadir, params.interpro.toString())
        }

        // Validate the interpro and Xref data files
        // Validate these first as we need the interpro files for validateAppData
        error = InterProScan.validateXrefFiles(_datadir, params.xRefsConfig, params.goterms, params.pathways)
        if (error) {
            log.error error
            exit 1
        }

        // Validate the selected member databases data files and dirs
        error = InterProScan.validateAppData(apps, _datadir, params.appsConfig)
        if (error) {
            log.error error
            exit 1
        }


    }

    if (need_data && params.download) {
        DOWNLOAD_INTERPRO(apps, params.appsConfig)
        _interproRelease = DOWNLOAD_INTERPRO.out.interproRelease.val
    } else if (need_data && !params.download) {
        // Check there is an InterPro dir
        (_interproDir, error) = InterProScan.validateInterproDir(_datadir)
        if (error) {  // no datadir/interpro dir found
            log.error error
            exit 1
        }

        if (params.interpro == "latest") {
            // Use the latest interpro release in datadir/interpro
            def _dirs = _interproDir.listFiles()
                        .findAll { it.isDirectory() }
                        .collect { it.name.toFloat } // what if the user put a dir in that can't be converted to a float
            if (_dirs.isEmpty()) {
                log.error ("No InterPro release directories were found in ${_datadir}/interpro"+
                           "Please ensure that the data dir is correctly populated or use --download")
               exit 1
            }
            _interproRelease = _dirs.max()
        } else {
            // Use the user specified InterPro release
            _interproRelease = params.interpro.toString()
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
        def _interproRelease = null
        def _interproDir = ""
        def _versionsMaps = new JsonSlurper().parse(new File(_versionsFile.toString())) // E.g. [6.0: [102.0, 103.0, 104.0]]

        if (!_versionsMaps.containsKey(_iprScanVersion)) {
            log.error "InterProScan version ${_iprScanVersion} is not listed in ${_versions}"
            exit 1
        }

        if (params.interpro == "latest") {
            def _compatibleVersions = _versionsMaps[_iprScanVersion]*.toFloat()
            def _interproDir = new File(_datadir.resolve("interpro").toString())

            if (!_interproDir.exists() || !_interproDir.isDirectory() {
                if (!params.download) {
                    log.error ("No 'interpro' directory was found in the data directory ${_datadir}\n" +
                               "Please ensure that the data dir is correctly populated or use --download")
                    exit 1
                }
                _interproRelease = _compatibleVersions.max()
                _downloadInterpro = true
                _downloadApps = true
                return
            }

            // Check which InterPro releases are already in the datadir, if any
            def _dirs = _interproDir.listFiles()
                        .findAll { it.isDirectory() }
                        .collect { it.name.toFloat }
                        .findAll { it in _compatibleVersions }

            if (_dirs.isEmpty()) {
                if (!params.download) {
                    log.error ("No compatible InterPro releases are stored in the data dir ${_datadir}.\n" +
                               "Compatible releases: ${_compatibleVersions.join(', ')}\n" +
                               "Please ensure that the data dir is correctly populated or use --download")
                    exit 1
                }
                _interproRelease = _compatibleVersions.max()
                _downloadInterpro = true
                _downloadApps = true
                return
            }

            def _latestDir = _dirs.max()
            // We already filtered the _dirs for compatible InterPro versions so _lastestDir > _compatibleVersions.max()
            // is impossible
            if (_latestDir < _compatibleVersions.max()) {
                if (params.download) {
                    _interproRelease = _compatibleVersions.max()
                    _downloadInterpro = true
                    _downloadApps = true
                } else {
                    _interproRelease = _dirs.max()
                    log.info("A newer compatible InterPro release (${_compatibleVersions.max()}) is available.\n" +
                             "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release")
                }
            } else {
                _interproRelease = _compatibleVersions.max()
            }
        } else {
            // Using the user selected InterPro release. But first check it is compatible with the IprScan version
            if (_versionsMaps[_iprScanVersion].contains(params.interpro.toString())) {
                _interproRelease = params.interpro.toString()

                // Check the interpro datadir exists
                def _interproDir = new File(_datadir.resolve("interpro").resolve(_interproRelease).toString())
                if (!_interproDir.exists() || !_interproDir.isDirectory() {
                    if (!params.download) {
                        log.error ("No 'interpro/${_interproRelease}' directory was found in the data directory ${_datadir}\n" +
                                   "Please ensure that the data dir is correctly populated or use --download")
                        exit 1
                    }
                    _downloadInterpro = true
                    _downloadApps = true
                    return
                }
            } else {
                log.error "Interpro version ${params.interpro} is not compatible with InterProScan version ${_iprScanVersion}"
                exit 1
            }
        }
        println "# InterPro: $_interproRelease"

        // [5] Check if need to download any member database data
        if (!_downloadApps) {
            missingApps = InterProScan.validateAppData(apps, _datadir, params.appsConfig, true)  // returns a Map of missing datafiles
            if (missingApps) {
                _downloadApps = true
            }
        }

        // [6] Download data if needed and enabled
        if (_downloadInterpro || _downloadApps) {
            DOWNLOAD_INTERPRO(apps, params.appsConfig, params.datadir, _interproRelease, _iprScanVersion)
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

def checkData() {
    /* This code is placed in a separate function so that if at any point download is set to true
    we exit and go straight to downloading data. */

}