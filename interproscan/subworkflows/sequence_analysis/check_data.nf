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

    def dataPaths = [:]
    boolean gene3d_funfam_processed = false

    applications.split(',').each { member ->

        if (member in ['antifam', 'ncbifam', 'smart', 'superfamily']) {
            dataPaths[member] = [hmm: "${dataDir}/${params.members."${member}".hmm}"]
        } else if (member == 'cdd') {
            dataPaths[member] = [
                library: "${dataDir}/${params.members."${member}".library}",
                postprocess_data: "${dataDir}/${params.members."${member}".postprocess.data}"
            ]
        } else if (member == 'gene3d' || member == 'funfam' && !gene3d_funfam_processed) {
            gene3d_funfam_processed = true
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."gene3d".hmm}",
                model2sf_map: "${dataDir}/${params.members."gene3d".postprocess.model2sf_map}",
                discontinuous_regs: "${dataDir}/${params.members."gene3d".postprocess.discontinuous_regs}",
                funfam_hmm: "${dataDir}/${params.members."funfam".hmm}"
            ]
        } else if (member == 'hamap') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                postprocess_models_dir: "${dataDir}/${params.members."${member}".postprocess.models_dir}"
            ]
        } else if (member == 'panther') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                postprocess_data_dir: "${dataDir}/${params.members."${member}".postprocess.data_dir}",
                postprocess_paint_annotations: "${dataDir}/${params.members."${member}".postprocess.paint_annotations}"
            ]
        } else if (member == 'pfam') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                postprocess_seed: "${dataDir}/${params.members."${member}".postprocess.seed}",
                postprocess_clan: "${dataDir}/${params.members."${member}".postprocess.clan}",
                postprocess_data: "${dataDir}/${params.members."${member}".postprocess.data}"
            ]
        } else if (member == 'pirsf') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                postprocess_data: "${dataDir}/${params.members."${member}".postprocess.data}"
            ]
        } else if (member == 'pirsr') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                postprocess_rules: "${dataDir}/${params.members."${member}".postprocess.rules}"
            ]
        } else if (member == 'prints') {
            dataPaths[member] = [
                hierarchy: "${dataDir}/${params.members."${member}".data.hierarchy}"
            ]
        } else if (member == "prosite_patterns") {
            dataPaths[member] = [
                data: "${dataDir}/${params.members."${member}".data}",
                evaluator: "${dataDir}/${params.members."${member}".evaluator}"
            ]
        } else if (member == "prosite_profiles") {
            dataPaths[member] = [
                data: "${dataDir}/${params.members."${member}".data}",
                skip_flagged_profiles: "${dataDir}/${params.members."${member}".skip_flagged_profiles}"
            ]
        } else if (member == 'sfld') {
            dataPaths[member] = [
                hmm: "${dataDir}/${params.members."${member}".hmm}",
                sites_annotation: "${dataDir}/${params.members."${member}".postprocess.sites_annotation}",
                hierarchy: "${dataDir}/${params.members."${member}".postprocess.hierarchy}"
            ]
        }
    }

    def missingFiles = []
    dataPaths.each { category, files ->
        files.each { key, filePath ->
            if (filePath.endsWith("Cdd_NCBI")) {
                /*
                There are multiple NCBI files
                Check at least some NCBI files are present
                */
                if (!file("${filePath}.pal").exists()) {
                    missingFiles << filePath
                }
            } else if (!file(filePath).exists()) {
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
