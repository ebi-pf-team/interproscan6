include { getMatchesApiUrl } from "../../modules/lookup"

workflow INIT_PIPELINE {
    // Validate pipeline input parameters

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

    emit:
    fasta            // str: path to input fasta file
    apps             // list: list of application to
    outdir           // str: path to output directory
    formats          // set<String>: output file formats
    signalpMode      // str: Models to be used with SignalP
    matchesApiUrl    // str|null: URL of Matches API to query
}
