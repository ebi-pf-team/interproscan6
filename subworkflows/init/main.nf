workflow INIT_PIPELINE {
    // Validate pipeline input parameters
    take:
    input
    applications
    apps_config
    runML
    datadir
    formats
    outdir
    outprefix
    no_matches_api
    matches_api_url
    interpro_version
    skip_intepro
    skip_applications
    goterms
    pathways
    workflow_manifest

    main:
    // Check the input
    fasta = InterProScan.resolveFile(input)
    if (!fasta) {
        log.error "No such file: ${input}"
        exit 1
    }

    // Applications validation
    (apps, error, warn) = InterProScan.validateApplications(applications, skip_applications, apps_config, runML)
    if (!apps) {
        log.error error
        exit 1
    } else if (warn) {
        log.warn warn
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
            exit 1
        }
    
        (datadir, error) = InterProScan.resolveDirectory(datadir, false, false)
        if (datadir == null) {
            log.error error
            exit 1
        }
    } else {
        datadir = null
    }
  
    version = InterProScan.validateInterProVersion(interpro_version)
    if (version == null) {
        log.error "--interpro <VERSION>: invalid format; expecting number or 'latest'"
        exit 1
    }

    (outdir, error) = InterProScan.resolveDirectory(outdir, false, false)
    if (!outdir) {
        log.error error
        exit 1
    }

    if (outprefix == null) {
        outprefix = "${outdir}/${fasta.split('/').last()}"
    } else if (outprefix.contains("/") || outprefix.contains(File.separator)) {
        log.error "--outprefix must not contain slashes or directory names. Use --outdir to control output location."
        exit 1
    } else {
        outprefix = "${outdir}/${outprefix}"
    }

    (matches_api_apps, local_only_apps, api_version, error) = Lookup.prepareLookup(
        apps,
        no_matches_api,
        matches_api_url,
        workflow_manifest
    )
    if (error) {
        log.warn error
    } else if (!no_matches_api && !local_only_apps.isEmpty()) {
        log.warn "The following analyses are not available in the Matches API: " +
                local_only_apps.join(", ") +
                ". They will be executed locally."
    }

    emit:
    fasta            // str: path to input fasta file
    local_only_apps  // list: list of application to that are not in the matches API
    matches_api_apps // list: list of applications that are in the matches API
    api_version      // str: version of the matches API
    datadir          // str: path to data directory, or null if not needed
    outprefix        // str: base path for output files
    formats          // set<String>: output file formats
    version          // str: InterPro version (or "latest")
}
