import groovy.json.JsonOutput


workflow CHECK_DATA {
    take:
    applications
    data_dir

    main:
    // check the data dir before check the paths individually
    def dDir = new File(data_dir)
    if (dDir.exists() && dDir.isDirectory()) {
        dataDir = data_dir.endsWith('/') ? data_dir[0..-2] : data_dir
    } else {
        log.error "Could not find data directory at '${data_dir}'.\nPlease check the value of --datadir"
        exit 5
    }

    def missingFiles = []
    boolean gene3d_funfam_processed = false

    user_applications = applications.toLowerCase()
    user_applications.split(',').each { member ->

        if (member in ['antifam', 'ncbifam', 'smart', 'superfamily']) {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'cdd') {
            filePath = "${dataDir}/${params.members."${member}".postprocess.data}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".library}"
            if (!file("${filePath}.pal").exists()) {
                missingFiles << filePath
            }
        } else if (member == 'gene3d' || member == 'funfam' && !gene3d_funfam_processed) {
            gene3d_funfam_processed = true
            filePath = "${dataDir}/${params.members."gene3d".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."gene3d".postprocess.model2sf_map}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."gene3d".postprocess.discontinuous_regs}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            if (user_applications.contains('funfam')) {
                filePath = "${dataDir}/${params.members."funfam".hmm}"
                if (!file(filePath).exists()) {
                    missingFiles << filePath
                }
            }
        } else if (member == 'hamap') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.models_dir}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'panther') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.data_dir}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.paint_annotations}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'pfam') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.seed}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.clan}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.data}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'pirsf') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.data}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'pirsr') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.rules}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'prints') {
            filePath = "${dataDir}/${params.members."${member}".data.hierarchy}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == "prosite_patterns") {
            filePath = "${dataDir}/${params.members."${member}".data}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".evaluator}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == "prosite_profiles") {
            filePath = "${dataDir}/${params.members."${member}".data}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".skip_flagged_profiles}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        } else if (member == 'sfld') {
            filePath = "${dataDir}/${params.members."${member}".hmm}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.sites_annotation}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
            filePath = "${dataDir}/${params.members."${member}".postprocess.hierarchy}"
            if (!file(filePath).exists()) {
                missingFiles << filePath
            }
        }
    }

    if (missingFiles) {
        log.error "Could not find all necessary data files in '${dataDir}/'\nMissing files:\n${missingFiles.join('\n')}"
        exit 5
    }

    emit:
    dataDir
}
