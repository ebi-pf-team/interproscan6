import groovy.json.JsonSlurper

include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../../modules/download"
include { FIND_MISSING_DATA             } from "../../../modules/download"
include { VALIDATE_DATA                 } from "../../../modules/download"


workflow PREPARE_DATABASES {
    take:
    applications
    apps_config
    data_dir
    interpro_version
    iprscan_version
    add_goterms
    add_pathways

    main:
    iprscan_major_minor = iprscan_version.split("\\.")[0..1].join(".")
    ch_ready = Channel.empty()

    if (data_dir != null) {
        // If data_dir is null, we only run analyses that do not depend on data files (e.g. coils)
        versions = InterProScan.fetchCompatibleVersions(iprscan_major_minor)
        if (versions == null) {
            log.warn """InterProScan could not retrieve compatibility information \
for InterPro data versions from EMBL-EBI FTP."""

            if (interpro_version == "latest") {
                highest_version = InterProScan.findLocalHighestVersionDir("${data_dir}/interpro")
                if (!highest_version) {
                    log.error "No version of InterPro found in ${data_dir}/interpro"
                    exit 1
                }

                interpro_version = highest_version
                log.warn """InterProScan is using the highest locally available InterPro \
data version (${interpro_version}), but compatibility with this version of InterProScan \
cannot be verified."""
            }
        } else if (interpro_version == "latest") {
            interpro_version = versions[-1]
        } else if (!versions.contains(interpro_version)) {
            error = "InterProScan ${iprscan_version} is not compatible "
            error += "with InterPro ${interpro_version} data.\n"
            error += "Compatible versions are: ${versions.join(', ')}."
            log.error error
            exit 1
        }

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
            data_dir
        )

        ch_ready = ch_ready.mix(DOWNLOAD_DATABASE.out)
        ch_ready = ch_ready.collect(flat: false)
    }

    VALIDATE_DATA(ch_ready.ifEmpty { [] })

    emit:
    versions = VALIDATE_DATA.out                // map: [ dbname: [version: <version>, path: <datapath>] ]
    iprscan_major_minor
}
