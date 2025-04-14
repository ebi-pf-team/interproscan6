
include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../modules/download"
include { PREPARE_DOWNLOADS             } from "../../modules/download"
include { GET_MEMBER_RELEASES           } from "../../modules/download"


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
           log.error "No version of InterPro found in ${datadir}/interpro"
           exit 1
        }

        interpro_version = highest_version
        log.warn """Without the '--download' option enabled, InterProScan uses \
the highest locally available version of InterPro data, but cannot \
verify compatibility with InterProScan. \
To ensure you're using the latest compatible data, re-run with the --download option."""
    }

    // Most members have a single dir, but gene3d and funfam are collated under cath for example
    app_dirs = apps_config
        .findAll { k, v -> v.has_data == true }
        .collectEntries { key, value ->
            def dir = value.get("dir", "")
            def parts = dir.split('/')
            def first = parts[0]
            def second = parts.size() > 1 ? parts[1..-1].join('/') : ""
            return [key, [first, second]]
    }

    // Check if we have the interpro database release file
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
        if (download) {
            DOWNLOAD_INTERPRO(
                ["interpro", interpro_version],
                iprscan_major_minor,
                datadir
            )
            // ready_ch = DOWNLOAD_INTERPRO.out.collect()
            // memberDbReleases = InterProScan.getMemberDbReleases(path, ready_ch)
            PREPARE_DOWNLOADS(
                DOWNLOAD_INTERPRO.out,
                path,
                applications,
                app_dirs,
                datadir
            )
        } else {
            log.error """No database release file found in ${datadir}/interpro/${interpro_version}
Use the '--download' option to automatically download InterPro release data."""
            exit 1
        }
    }

    ch_appls = PREPARE_DOWNLOADS.out.flatMap()

    // Create an empty Channel to be used as default when no databases need download
    empty_db_channel = Channel.of([]).collect()

    if (download) {
        ch_ready = DOWNLOAD_DATABASE(
            ch_appls,
            iprscan_major_minor,
            datadir
        ).collect()
    } else {
        // Check if we have apps that need download
        ch_appls_check = ch_appls.collect()
        ch_appls_check.subscribe { apps ->
            if (apps.size() > 0) {
                def msg = apps.collect { app -> "Database: ${app[0]} - Release: ${app[1]}" }.join("\n")
                log.error """Data is missing in ${datadir} for the following applications:
${msg}
Use the '--download' option to automatically download InterPro release data."""
                exit 1
            }
        }
    }

    // Force wait on the databases.json path whether we have apps or not
    memberDbReleases = GET_MEMBER_RELEASES(
        path,
        ch_appls.collect().ifEmpty { empty_db_channel }
    )
    interproscanVersion = iprscan_major_minor

    emit:
    memberDbReleases       // map: [db name (lowercase): release]
    interproscanVersion
}
