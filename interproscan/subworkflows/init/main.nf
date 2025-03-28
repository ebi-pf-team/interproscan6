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
            if (!params.datadir) {
                log.error "--datadir <DATA-DIR> is mandatory when using --download"
                exit 1
            }
            DOWNLOAD_INTERPRO(apps, params.appsConfig)
            _interproRelease = DOWNLOAD_INTERPRO.out.interproRelease.val
        } else {
            // Check if there is a data directory
            // If --datadir is called and no path is given it converts to a boolean
            def dirPath = params.datadir instanceof Boolean ? null : params.datadir
            (datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
            if (!datadir) {
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
                (noConnWarning, latestReleaseWarning, compatError) = InterPro.checkCompatibility(_iprScanVersion, _interproRelease)
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
        println "# InterPro: $_interproRelease"  // adds to the pipeline intro lines in main.nf

        /* TODO: Update InterProScan.validate<TYPE>Data/Files
        to include the release numbers. The applciations.conf has 
        already been updated. */

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