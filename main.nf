include { INIT_PIPELINE } from './subworkflows/init/pipeline'
include { INTERPROSCAN  } from './workflows/interproscan'

params {
    help: Boolean
    globus: Boolean
    noMatchesApi: Boolean
    nucleic: Boolean
    goterms: Boolean
    pathways: Boolean
    skipReprLocations: Boolean
    skipInterproVersionCheck: Boolean
    runMl: Boolean
    useGpu: Boolean
    batchSize: Integer
    subBatchSize: Integer
    cpus: Integer
    matchesApiChunkSize: Integer
    matchesApiMaxRetries: Integer
    maxWorkers: Integer
    localCopyDir: Path
    input: Path
    datadir: Path
    outdir: Path
    applications: String
    formats: String
    outprefix: String
    skipApplications: String
    interpro: String
    matchesApiUrl: String
}

workflow {
    println "# ${workflow.manifest.name} ${workflow.manifest.version}"
    println "# ${workflow.manifest.description}\n"

    if (params.help) {
        uk.ac.ebi.interpro.InterProScan.printHelp(params.appsConfig)
        exit 0
    }

    uk.ac.ebi.interpro.InterProScan.validateParams(params, log)

    INIT_PIPELINE(
        params.input,
        params.applications,
        params.appsConfig,
        params.runMl,
        params.useGpu,
        params.datadir,
        params.formats,
        params.outdir,
        params.outprefix,
        params.interpro,
        params.skipReprLocations,
        params.skipApplications,
        params.goterms,
        params.pathways,
        workflow.manifest
    )
    fasta_file           = channel.fromPath(INIT_PIPELINE.out.fasta.val)
    apps                 = INIT_PIPELINE.out.apps.val
    apps_config          = INIT_PIPELINE.out.apps_config.val
    data_dir             = INIT_PIPELINE.out.datadir.val
    outprefix            = INIT_PIPELINE.out.outprefix.val
    formats              = INIT_PIPELINE.out.formats.val
    interpro_version     = INIT_PIPELINE.out.version.val

    INTERPROSCAN(
        fasta_file,
        apps,
        apps_config,
        data_dir,
        outprefix,
        formats,
        interpro_version,
        workflow.manifest.version,
        workflow.manifest.name,
        params.noMatchesApi,
        params.matchesApiUrl,
        params.matchesApiChunkSize,
        params.matchesApiMaxRetries,
        params.batchSize,
        params.subBatchSize,
        params.nucleic,
        params.skipReprLocations,
        params.goterms,
        params.pathways,
        params.globus,
        !params.skipInterproVersionCheck
    )

    workflow.onComplete = {
        if (workflow.success) {
            println "\nInterProScan completed successfully"
            println "Results available at: ${outprefix}.*"
            if (workflow.duration.toSeconds() <= 60) {
                def formatter = java.time.format.DateTimeFormatter.ofPattern("dd-MMM-yyyy HH:mm:ss");
                println "Completed at        : ${workflow.complete.format(formatter)}"
                println "Duration            : ${workflow.duration}"
            }
        }        
    }
}
