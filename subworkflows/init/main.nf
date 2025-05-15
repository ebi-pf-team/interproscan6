workflow INIT_PIPELINE {
    // Validate pipeline input parameters
    take:
    input
    applications
    apps_config
    download
    offline
    datadir
    formats
    outdir
    matches_api_url
    interpro_version
    skip_intepro
    goterms
    pathways

    main:
    // Check the input
    fasta = InterProScan.resolveFile(input)
    if (!fasta) {
        log.error "No such file: ${input}"
        exit 1
    }

    // Applications validation
    (apps, error) = InterProScan.validateApplications(applications, apps_config)
    if (!apps) {
        log.error error
        exit 1
    }

    if (download && offline) {
        log.error "--download and --offline are mutually exclusive"
        exit 1
    }

    if (skip_intepro && (goterms || pathways)) {
        log.error "--skip_intepro is mutually exclusive with --goterms and --pathways"
        exit 1
    }

    // Check valid output file formats were provided
    (formats, error) = InterProScan.validateFormats(formats)
    if (error) {
        log.error error
        exit 1
    }

    apps_with_data = InterProScan.getAppsWithData(apps, apps_config)
    if (apps_with_data.size() > 0) {
        if (datadir == null) {
            log.error "'--datadir <DATA-DIR>' is required for the selected applications."
            edit 1
        }
    
        (datadir, error) = InterProScan.resolveDirectory(datadir, false, true)
        if (datadir == null) {
            log.error error
            exit 1
        }
    } else {
        datadir = null
    }
  
    version = InterProScan.validateInterProVersion(interpro_version)
    if (version == null) {
        log.error "--interpro <VERSION>: invalid format; expecting number of 'latest'"
        exit 1
    }

    (outdir, error) = InterProScan.resolveDirectory(outdir, false, true)
    if (!outdir) {
        log.error error
        exit 1
    }

    if (!offline) {
        invalidApps = apps.findAll { app ->
            ["signalp_euk", "signalp_prok", "deeptmhmm"].contains(app)
        }

        if (invalidApps) {
            log.error "Pre-calculated results for DeepTMHMM, SignalP_Euk, and SignalP_Prok are not yet available in the Matches API. To ensure these analyses run locally and produce results, please add the '--offline' flag when invoking the pipeline."
            exit 1
        }
    }

    emit:
    fasta            // str: path to input fasta file
    apps             // list: list of application to
    datadir          // str: path to data directory, or null if not needed
    outdir           // str: path to output directory
    formats          // set<String>: output file formats
    version          // str: InterPro version (or "latest")
}
