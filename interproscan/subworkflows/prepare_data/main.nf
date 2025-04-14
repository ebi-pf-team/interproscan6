import groovy.json.JsonSlurper

include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../modules/download"
include { VALIDATE_APP_DATA             } from "../../modules/download"
include { GET_DB_RELEASES           } from "../../modules/download"


workflow PREPARE_DATA {
    take:
    applications
    apps_config
    xref_config
    data_dir
    interpro_version
    iprscan_version
    download
    add_goterms
    add_pathways

    main:
    def iprscan_major_minor = iprscan_version.split("\\.")[0..1].join(".")

    /* Do not skip this step, even when only running members without data
    because we need the interpro version to check against the match lookup service
    API to make sure they are compatible. */
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
        def highest_version = InterProScan.findLocalHighestVersionDir("${data_dir}/interpro")
        if (!highest_version) {
           log.error "No version of InterPro found in ${data_dir}/interpro"
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
    db_json_path = "${data_dir}/${xref_config.dir}/${interpro_version}/databases.json"
    if (InterProScan.resolveFile(db_json_path)) {
        VALIDATE_APP_DATA(
            true,
            db_json_path,
            applications,
            app_dirs,
            data_dir
        )
    } else {
        if (download) {
            DOWNLOAD_INTERPRO(
                ["interpro", interpro_version],
                iprscan_major_minor,
                data_dir
            )

            VALIDATE_APP_DATA(
                DOWNLOAD_INTERPRO.out,
                db_json_path,
                applications,
                app_dirs,
                data_dir
            )
        } else {
            log.error """No database release file found in ${data_dir}/interpro/${interpro_version}
Use the '--download' option to automatically download InterPro release data."""
            exit 1
        }
    }

    // TODO: Fix this warning:: WARN: Input tuple does not match tuple declaration in process `PREPARE_DATA:DOWNLOAD_DATABASE` -- offending value: DataflowQueue(queue=[]
    ch_appls = VALIDATE_APP_DATA.out.map { jsonFile ->
        def jsonSlurper = new JsonSlurper()
        def data = jsonSlurper.parse(new File(jsonFile.toString()))
        return data
    }.flatMap { items ->
        if (items.isEmpty()) {
            return Channel.empty()
        } else {
            return Channel.from(items.toList())
        }
    }

    // Create an empty Channel to be used as default when no databases have data that needs downloading
    empty_db_channel = Channel.of([]).collect()

    if (download) {
        ch_ready = DOWNLOAD_DATABASE(
            ch_appls,
            iprscan_major_minor,
            data_dir
        ).collect()
    } else {
        // Check if we have apps that need download
        ch_appls_check = ch_appls.collect()
        ch_appls_check.subscribe { apps ->
            if (apps.size() > 0) {
                def msg = apps.collect { app -> "Database: ${app[0]} - Release: ${app[1]}" }.join("\n")
                log.error """Data is missing in ${data_dir} for the following applications:
${msg}
Use the '--download' option to automatically download InterPro release data."""
                exit 1
            }
        }
    }

    // Force wait on the databases.json path whether we have apps or not
    path = "${data_dir}/${xref_config.dir}/${interpro_version}"
    memberDbReleasesPath = GET_DB_RELEASES(
        db_json_path,
        path,
        xref_config,
        add_goterms,
        add_pathways,
        ch_appls.collect().ifEmpty { empty_db_channel }
    )

    // Wait for the channel to resolve and assign the value
    memberDbReleasesPath.map { dbFilePath ->
        def jsonSlurper = new JsonSlurper()
        return jsonSlurper.parse(new File(dbFilePath.toString()))
    }.set { memberDbReleases }

    interproscanVersion = iprscan_major_minor

    emit:
    memberDbReleases       // map: [db name (lowercase): release]
    interproscanVersion
}
