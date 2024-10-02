workflow INIT_PIPELINE {
    main:
    // Params validation
    InterProScan.validateParams(params, log)
    fasta = InterProScan.resolveFile(params.input)
    if (!fasta) {
        log.error "No such file: ${params.input}"
        exit 1
    }
    (datadir, error) = InterProScan.resolveDirectory(params.datadir, true, false)
    if (!datadir) {
        log.error error
        exit 1
    }
    (outdir, error) = InterProScan.resolveDirectory(params.outdir, false, true)
    if (!outdir) {
        log.error error
        exit 1
    }
    
    // Applications validation
    (apps, error) = InterProScan.validateApplications(params.applications, params.appsConfig)
    if (!apps) {
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

    emit:
    fasta      // str: path to input fasta file
    datadir    // str: path to data directory
    apps       // list: list of application to
    outdir     // str: path to output directory
}