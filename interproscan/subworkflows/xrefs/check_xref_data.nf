import groovy.json.JsonOutput

workflow CHECK_XREF_DATA {
    take:
    dataDir

    main:
    def missingXrefsFiles = []

    def checkFilePath = { path ->
        if (!file(path).exists()) {
            missingXrefsFiles << path
        }
    }

    checkFilePath("${dataDir}/${params.xrefs.entries}")
    checkFilePath("${dataDir}/${params.xrefs.goterms}.json")
    checkFilePath("${dataDir}/${params.xrefs.goterms}.ipr.json")
    checkFilePath("${dataDir}/${params.xrefs.pathways}.json")
    checkFilePath("${dataDir}/${params.xrefs.pathways}.ipr.json")

    if (missingXrefsFiles) {
        log.error "Could not find all necessary XREF data files in '${dataDir}/'\nMissing XREF files:\n${missingXrefsFiles.join('\n')}"
        exit 5
    }
}
