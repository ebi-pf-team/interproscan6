import groovy.json.JsonSlurper

include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../../modules/download"
include { FIND_MISSING_DATA             } from "../../../modules/download"
include { VALIDATE_DATA                 } from "../../../modules/download"


workflow PREPARE_DATABASES {
    take:
    local_only_apps
    matches_api_apps
    apps_config
    data_dir
    interpro_version
    iprscan_version
    no_matches_api
    add_goterms
    add_pathways
    use_globus

    main:
    applications = local_only_apps + matches_api_apps
    iprscan_major_minor = iprscan_version.split("\\.")[0..1].join(".")
    ch_ready = Channel.empty()

    if (data_dir == null && no_matches_api) {
        /*
        If data_dir is not specified, we only run analyses that do not depend on data files (e.g. coils).
        We also don't need the InterPro version as we won't be using the Matches API either.
        But we still need to create a dummy channel for the output of VALIDATE_DATA to not be empty.
        We don't need the InterPro data dir
        */
        ch_ready = Channel.of(["default", "1.0", null])
    } else if (data_dir == null && !no_matches_api) {
        /*
        We don't have a datadir so don't check to download data
        But we still want to use the Matches API so we need the InterPro version
        We don't need the InterPro data dir
        */
        interpro_version = getInterproVersion(interpro_version, iprscan_major_minor, use_globus)
        ch_ready = Channel.of(["interpro", interpro_version, null])
    } else {
        /*
        Check the InterPro release version is compatible with the IPRScan version
        and if data files are needed and missing, download them
        */
        interpro_version = getInterproVersion(interpro_version, iprscan_major_minor, use_globus)

        // Most members have a single dir, but CATH-Gene3D and CATH-FuNFam are collated under cath for example
        app_dirs = apps_config
            .findAll { k, v -> v.has_data == true }
            .collectEntries { key, value ->
                def dir = value.get("dir", "")
                def parts = dir.split('/')
                def first = parts[0]
                def second = parts.size() > 1 ? parts[1..-1].join('/') : ""
                return [key, [first, second]]
        }

        // Not all members need data, if none of the applications need data skip downloading
        db_json_path = "${data_dir}/interpro/${interpro_version}/databases.json"
        if (InterProScan.resolveFile(db_json_path)) {
            // JSON file of database metadata found
            FIND_MISSING_DATA(
                ["", "", ""],  // state dependency, can be anything
                db_json_path,
                applications,
                app_dirs,
                data_dir
            )

            ch_interpro = Channel.value(["interpro", interpro_version, "${data_dir}/interpro/${interpro_version}"])
        } else {
            // Not found: download the InterPro metadata archive
            DOWNLOAD_INTERPRO(
                ["interpro", "interpro", interpro_version, false, "${data_dir}/interpro/${interpro_version}"],
                iprscan_major_minor,
                use_globus,
                data_dir
            )

            ch_interpro = DOWNLOAD_INTERPRO.out

            FIND_MISSING_DATA(
                ch_interpro,
                db_json_path,
                applications,
                app_dirs,
                data_dir
            )
        }

        ch_ready = ch_ready.mix(ch_interpro)
        ch_ready = ch_ready.mix(FIND_MISSING_DATA.out.with_data.flatMap())
        ch_to_download = FIND_MISSING_DATA.out.without_data.flatMap()
    
        DOWNLOAD_DATABASE(
            ch_to_download,
            iprscan_major_minor,
            use_globus,
            data_dir
        )

        ch_ready = ch_ready.mix(DOWNLOAD_DATABASE.out)
    }

    ch_ready = ch_ready.collect(flat: false)
    VALIDATE_DATA(ch_ready)

    emit:
    versions = VALIDATE_DATA.out                // map: [ dbname: [version: <version>, path: <datapath>] ]
    iprscan_major_minor
}

def getInterproVersion(String interpro_version, String iprscan_major_minor, boolean use_globus) {
    versions = InterProScan.fetchCompatibleVersions(iprscan_major_minor, use_globus)
    if (versions == null) {
        if (use_globus) {
            log.warn "InterProScan could not retrieve compatibility information" +
                        " for InterPro data versions."
        } else {
            log.warn "InterProScan could not retrieve compatibility information" +
                        " for InterPro data versions from the EMBL-EBI FTP.\n" +
                        "Try using the --globus option to use the Globus mirror instead."
        }

        if (interpro_version == "latest") {
            highest_version = InterProScan.findLocalHighestVersionDir("${data_dir}/interpro")
            if (!highest_version) {
                log.error "No version of InterPro found in ${data_dir}/interpro"
                exit 1
            }

            interpro_version = highest_version
            log.warn "InterProScan is using the highest locally available InterPro" +
                        " data version (${interpro_version}), but compatibility with" +
                        " this version of InterProScan cannot be verified."
        }
    } else if (interpro_version == "latest") {
        interpro_version = versions[-1]
    } else if (!versions.contains(interpro_version)) {
        error = "InterProScan ${iprscan_version} is not compatible with InterPro " +
                " ${interpro_version} data.\n" +
                "Compatible versions are: ${versions.join(', ')}."
        log.error error
        exit 1
    }

    return interpro_version
}
