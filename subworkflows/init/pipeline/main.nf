workflow INIT_PIPELINE {
    // Validate pipeline input parameters
    take:
    input
    applications
    applications_config
    run_ml
    use_gpu
    datadir
    formats
    outdir
    outprefix
    seqdbdir
    interpro_version
    skip_intepro
    skip_applications
    goterms
    pathways
    _workflow_manifest

    main:
    // Check the input
    fasta = file(input)
    if (!fasta.isFile()) {
        log.error "No such file: ${input}"
        exit 1
    }

    // Applications validation
    (apps, error, warn) = uk.ac.ebi.interpro.InterProScan.validateApplications(applications, skip_applications, applications_config, run_ml)
    if (!apps) {
        log.error error
        exit 1
    } else if (warn) {
        log.warn warn
    }

    // Enable gpu acceleration if requested
    (apps_config, warn) = uk.ac.ebi.interpro.InterProScan.enableGpuAcceleration(use_gpu, apps, applications_config)
    if (warn) {
        log.warn warn
    }

    if (skip_intepro && (goterms || pathways)) {
        log.error "--skip_intepro is mutually exclusive with --goterms and --pathways"
        exit 1
    }

    // Check valid output file formats were provided
    (formats, error) = uk.ac.ebi.interpro.InterProScan.validateFormats(formats)
    if (error) {
        log.error error
        exit 1
    }

    apps_with_data = uk.ac.ebi.interpro.InterProScan.getAppsWithData(apps, apps_config)
    if (apps_with_data.size() > 0) {
        if (datadir == null) {
            log.error "'--datadir <DATA-DIR>' is required for the selected applications."
            exit 1
        }

        datadir = file(datadir)
        if (datadir.isFile()) {
            log.error "'--datadir <DATA-DIR>' is required and cannot be an existing file."
            exit 1
        } else if (!datadir.isDirectory()) {
            datadir.mkdirs()
        }            
    } else {
        datadir = null
    }
  
    version = uk.ac.ebi.interpro.InterProScan.validateInterProVersion(interpro_version)
    if (version == null) {
        log.error "--interpro <VERSION>: invalid format; expecting number or 'latest'"
        exit 1
    }

    if (outdir == null) {
        log.error "'--outdir <OUTPUT-DIR>' is required."
        exit 1
    }

    outdir = file(outdir)
    if (outdir.isFile()) {
        log.error "'--outdir <OUTPUT-DIR>' is required and cannot be an existing file."
        exit 1
    } else if (!outdir.isDirectory()) {
        outdir.mkdirs()
    }

    if (outprefix == null) {
        outprefix = outdir.resolve(fasta.name)
    } else if (outprefix.contains("/") || outprefix.contains(File.separator)) {
        log.error "--outprefix must not contain slashes or directory names. Use --outdir to control output location."
        exit 1
    } else {
        outprefix = outdir.resolve(outprefix)
    }

    if (seqdbdir != null) {
        seqdbdir = file(seqdbdir)
        if (seqdbdir.isFile()) {
            log.error "'--seqdbdir <SEQDBDIR>' cannot be an existing file."
            exit 1
        } else if (!seqdbdir.isDirectory()) {
            seqdbdir.mkdirs()
        }
    }

    emit:
    fasta            // path: path to input fasta file
    apps             // list: selected applications
    apps_config      // map: updated applications configuration
    datadir          // path: path to data directory, or null if not needed
    outprefix        // path: base path for output files
    formats          // set<String>: output file formats
    version          // str: InterPro version (or "latest")
}
