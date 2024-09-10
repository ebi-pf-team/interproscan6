import groovy.json.JsonOutput

workflow CHECK_XREF_DATA {
    take:
    dataDir
    goterms
    pathways

    main:
    def missingXrefsFiles = []

    def checkFilePath = { path ->
        if (!file(path).exists()) {
            missingXrefsFiles << path
        }
    }

    checkFilePath("${dataDir}/${params.xrefs.entries}")
    if (goterms) {
        checkFilePath("${dataDir}/${params.xrefs.goterms}.json")
        checkFilePath("${dataDir}/${params.xrefs.goterms}.ipr.json")
    }
    if (pathways) {
        checkFilePath("${dataDir}/${params.xrefs.pathways}.json")
        checkFilePath("${dataDir}/${params.xrefs.pathways}.ipr.json")
    }
    
    xrefDataDir = params.xrefs.entries.substring(0, params.xrefs.entries.lastIndexOf('/'))

    emit:
        missingXrefs = missingXrefsFiles.join('\n')
        xrefDataDir
}
