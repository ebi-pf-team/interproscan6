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

    // Helper function to check file existence and add to missingFiles list
    def checkFilePath = { path ->
        if (!file(path).exists()) {
            missingFiles << path
        }
    }

    user_applications = applications.toLowerCase()
    user_applications.split(',').each { member ->
        if (member in ['antifam', 'ncbifam', 'smart', 'superfamily']) {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
        } else if (member == 'cdd') {
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.data}")
            checkFilePath("${dataDir}/${params.members."${member}".library}.pal")
        } else if (member == 'gene3d' || member == 'funfam' && !gene3d_funfam_processed) {
            gene3d_funfam_processed = true
            checkFilePath("${dataDir}/${params.members."gene3d".hmm}")
            checkFilePath("${dataDir}/${params.members."gene3d".postprocess.model2sf_map}")
            checkFilePath("${dataDir}/${params.members."gene3d".postprocess.discontinuous_regs}")
            if (user_applications.contains('funfam')) {
                checkFilePath("${dataDir}/${params.members."funfam".hmm}")
            }
        } else if (member == 'hamap') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.models_dir}")
        } else if (member == 'panther') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.data_dir}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.paint_annotations}")
        } else if (member == 'pfam') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.seed}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.clan}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.data}")
        } else if (member == 'pirsf') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.data}")
        } else if (member == 'pirsr') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.rules}")
        } else if (member == 'prints') {
            checkFilePath("${dataDir}/${params.members."${member}".data.hierarchy}")
        } else if (member == "prosite_patterns") {
            checkFilePath("${dataDir}/${params.members."${member}".data}")
            checkFilePath("${dataDir}/${params.members."${member}".evaluator}")
        } else if (member == "prosite_profiles") {
            checkFilePath("${dataDir}/${params.members."${member}".data}")
            checkFilePath("${dataDir}/${params.members."${member}".skip_flagged_profiles}")
        } else if (member == 'sfld') {
            checkFilePath("${dataDir}/${params.members."${member}".hmm}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.sites_annotation}")
            checkFilePath("${dataDir}/${params.members."${member}".postprocess.hierarchy}")
        }
    }

    if (missingFiles) {
        log.error "Could not find all necessary data files in '${dataDir}/'\nMissing files:\n${missingFiles.join('\n')}"
        exit 5
    }

    emit:
    dataDir
}
