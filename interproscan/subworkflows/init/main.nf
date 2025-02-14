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

    // Validate the data dir and application data files if needed by any members
    // e.g. mobidblite and coils do no need additional data files
    need_data = apps.any { String appName ->
        params.appsConfig.get(appName)?.has_data || InterProScan.LICENSED_SOFTWARE.contains(appName)
    }
    if (need_data) {
        // If --datadir is called and no path is given it converts to a boolean
        def dirPath = params.datadir instanceof Boolean ? null : params.datadir
        (datadir, error) = InterProScan.resolveDirectory(dirPath, true, false)
        if (!datadir) {
            log.error error
            exit 1
        }
        error = InterProScan.validateAppData(apps, datadir, params.appsConfig)
        if (error) {
            log.error error
            exit 1
        }
        error = InterProScan.validateXrefFiles(datadir, params.xRefsConfig, params.goterms, params.pathways)
        if (error) {
            log.error error
            exit 1
        }
    } else {
        datadir = ""
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

    // Create the path to the database. The database must be built in a process running
    // within the ips6 container, thus the database is built in populate_database
    // String dbPath = "${workflow.workDir}/ips6.seq.db"

    emit:
    fasta        // str: path to input fasta file
    datadir      // str: path to data directory
    apps         // list: list of application to
    outdir       // str: path to output directory
    formats      // set<String>: output file formats
    signalpMode  // str: Models to be used with SignalP
}