// codenarc-disable ModuleIncludedTwiceRule
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.InterProScan

include { DOWNLOAD as DOWNLOAD_INTERPRO } from "../../../modules/download"
include { DOWNLOAD as DOWNLOAD_DATABASE } from "../../../modules/download"
include { FIND_DATABASES                } from "../../../modules/download"


workflow INIT_DATABASES {
    take:
    applications            // list[str] of applications
    appl_configs            // map of application config
    data_dir                // path
    interpro_version        // str, version of InterPro to run
    iprscan_version         // str, version of InterProScan
    use_matches_api         // bool, whether to use the Matches API
    use_globus              // bool, whether to use Globus instead of the EMBL-EBI FTP
    enforce_compatibility   // bool, whether we check that InterProScan and InterPro are compatible

    main:
    iprscan_maj_min_version = extractMajorMinorVersion(iprscan_version)

    if (interpro_version == "latest" || enforce_compatibility) {
        // if (interpro_version == "latest" && enforce_compatibility) {
        //     log.warn "--skip-interpro-version-check and --interpro 'latest' are mutually exclusive."
        //     enforce_compatibility = true
        // }

        def versions = InterProScan.fetchCompatibleVersions(iprscan_maj_min_version, use_globus)
        if (versions == null) {
            def message = null
            if (use_globus) {
                message = "InterProScan could not retrieve compatibility information for InterPro data versions. Try disabling the compatibility check with --skip-interpro-version-check."
            } else {
                message = "InterProScan could not retrieve compatibility information for InterPro data versions from the EMBL-EBI FTP. Try using the Globus mirror with --globus or disabling the compatibility check with --skip-interpro-version-check."
            }

            log.error message
            exit 1
        } else if (interpro_version == "latest") {
            interpro_version = versions[-1]
        } else if (!versions.contains(interpro_version)) {
            log.error "InterProScan ${iprscan_version} is not compatible with InterPro ${interpro_version} data. Compatible versions are: ${versions.join(', ')}."  // codenarc-disable-line JoinMismatchRule, JoinDuplicateRule
            exit 1
        }
    }

    ch_ready = Channel.empty()
    if (data_dir != null) {
        // At least one applications requires data files
        
        // Handle applications with a common directory (cath -> CATH-Gene3 / CATH-FunFam, prosite -> PROSITE Patterns / PROSITE Profiles)
        appl_dirs = appl_configs
            .findAll { k, v -> v.has_data == true }
            .collectEntries { key, value ->
                def dir = value.get("dir", "")
                def parts = dir.split('/')
                def first = parts[0]
                def second = parts.size() > 1 ? parts[1..-1].join('/') : ""  // codenarc-disable-line JoinMismatchRule, JoinDuplicateRule
                return [key, [first, second]]
        }

        // Find the JSON file for database metadata
        def interpro_dir = data_dir.resolve("interpro/${interpro_version}")
        def databases_json = interpro_dir.resolve("databases.json")
        if (databases_json.isFile()) {
            // InterPro data already downloaded (at least databases.json)
            ch_interpro = Channel.value(["interpro", interpro_dir])

            FIND_DATABASES(
                ["", ""],  // state dependency, can be anything
                databases_json,
                applications,
                appl_dirs,
                data_dir
            )
        } else {
            // Download the InterPro metadata archive
            ch_interpro = DOWNLOAD_INTERPRO(
                ["interpro", "interpro", interpro_version, false, interpro_dir],
                iprscan_maj_min_version,
                use_globus,
                data_dir
            )

            FIND_DATABASES(
                ch_interpro,
                db_json_path,
                applications,
                app_dirs,
                data_dir
            )
        }

        ch_ready = ch_ready.mix(ch_interpro)
        ch_ready = ch_ready.mix(FIND_DATABASES.out.ready.flatMap())
        ch_missing = FIND_DATABASES.out.missing.flatMap()
    
        DOWNLOAD_DATABASE(
            ch_missing,
            iprscan_maj_min_version,
            use_globus,
            data_dir
        )

        ch_ready = ch_ready.mix(DOWNLOAD_DATABASE.out)
    } else {
        ch_ready = Channel.of(["interpro", null])
    }

    emit:
    directories = ch_ready
    interpro_version
}


def extractMajorMinorVersion(String version) {
    return version.split("\\.")[0..1].join(".")  // codenarc-disable-line JoinMismatchRule, JoinDuplicateRule
}