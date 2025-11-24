nextflow.enable.dsl=2
import java.time.format.DateTimeFormatter
import uk.ac.ebi.interpro.InterProScan

include { INIT_PIPELINE      } from "./subworkflows/init"
include { INTERPROSCAN       } from "./workflows/interproscan.nf"

workflow {
    println "# ${workflow.manifest.name} ${workflow.manifest.version}"
    println "# ${workflow.manifest.description}\n"

    if (params.keySet().any { it.equalsIgnoreCase("help") }) {
        InterProScan.printHelp(params.appsConfig)
        exit 0
    }

    InterProScan.validateParams(params, log)

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
        params.skipInterpro,
        params.skipApplications,
        params.goterms,
        params.pathways,
        workflow.manifest
    )
    fasta_file           = Channel.fromPath(INIT_PIPELINE.out.fasta.val)
    apps                 = INIT_PIPELINE.out.apps.val
    apps_config          = INIT_PIPELINE.out.apps_config.val
    data_dir             = INIT_PIPELINE.out.datadir.val
    outprefix            = INIT_PIPELINE.out.outprefix.val
    formats              = INIT_PIPELINE.out.formats.val
    interpro_version     = INIT_PIPELINE.out.version.val

    INTERPROSCAN(
        fasta_file,
        apps,
        params.noMatchesApi,
        params.matchesApiUrl,
        apps_config,
        data_dir,
        outprefix,
        formats,
        interpro_version,
        workflow.manifest.version,
        workflow.manifest.name
    )

    workflow.onComplete = {
        if (workflow.success) {
            println "\nInterProScan completed successfully"
            println "Results available at: ${outprefix}.*"
            if (workflow.duration.toSeconds() <= 60) {
                DateTimeFormatter formatter = DateTimeFormatter.ofPattern("dd-MMM-yyyy HH:mm:ss");
                println "Completed at        : ${workflow.complete.format(formatter)}"
                println "Duration            : ${workflow.duration}"
            }
        }        
    }
}
