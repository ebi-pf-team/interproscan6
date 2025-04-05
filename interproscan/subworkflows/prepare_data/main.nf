
include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../modules/download"
include { PREPARE_DOWNLOADS             } from "../../modules/download"


workflow PREPARE_DATA {
    take:
    applications
    apps_config
    interpro_version
    iprscan_version
    datadir
    download

    main:
    def iprscan_major_minor = iprscan_version.split("\\.")[0..1].join(".")

    // TODO: need to bypass this if running application without data

    if (download) {
        def versions = InterProScan.fetchCompatibleVersions(iprscan_major_minor)
        if (versions == null) {
            log.error "Failed to retrieve compatible InterPro versions from EMBL-EBI FTP."
            exit 1
        } else if (interpro_version == "latest") {
            interpro_version = versions[-1]
        } else if (!versions.contains(interpro_version)) {
            def error = "InterProScan ${iprscan_version} is not compatible "
            error += "with InterPro ${interpro_version} data.\n"
            error += "Compatible versions are: ${versions.join(', ')}."
            log.error error
        }
    } else if (interpro_version == "latest") {
        def highest_version = InterProScan.findLocalHighestVersionDir("${datadir}/interpro")
        if (!highest_version) {
           log.error "No version of InterPro found"
           exit 1
        }
        
        interpro_version = highest_version
        log.warn """without the '--download' option enabled, InterProScan uses \
the highest locally available version of InterPro data, but cannot \
verify compatibility with InterProScan. \
To ensure you're using the latest compatible data, re-run with the --download option."""
    }
    
    app_dirs = apps_config
        .findAll { k, v -> v.has_data == true }
        .collectEntries { key, value ->
            def dir = value.get("dir", "")
            def parts = dir.split('/')
            def first = parts[0]
            def second = parts.size() > 1 ? parts[1..-1].join('/') : ""
            return [key, [first, second]]
    }

    path = "${datadir}/interpro/${interpro_version}/databases.json"
    if (InterProScan.resolveFile(path)) {
        PREPARE_DOWNLOADS(
            true,
            path,
            applications,
            app_dirs,
            datadir
        )
    } else {
        DOWNLOAD_INTERPRO(
            ["interpro", interpro_version],
            iprscan_major_minor,
            datadir
        )

        PREPARE_DOWNLOADS(
            DOWNLOAD_INTERPRO.out,
            path,
            applications,
            app_dirs,
            datadir
        )
    }

    ch_appls = PREPARE_DOWNLOADS.out.flatMap()

    ch_ready = DOWNLOAD_DATABASE(
        ch_appls,
        iprscan_major_minor,
        datadir
    ).collect()

    ch_ready.view()

    // emit:
    // interproRelease       // str: InterPro release version number
}
