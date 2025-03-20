import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'local'

    input:
    tuple val(index), val(fasta)
    val applications
    val url
    val chunkSize
    val maxRetries

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), optional: true

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def calculatedMatches = [:]

    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")
    def noLookupFasta = new StringBuilder()

    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]
    def md5List = sequences.keySet().toList()
    def chunks = md5List.collate(chunkSize)

    String baseUrl = InterPro.sanitizeURL(url.toString())
    boolean success = true
    for (chunk in chunks) {
        String data = JsonOutput.toJson([md5: chunk])
        def response = InterPro.httpRequest("${baseUrl}/matches", data, maxRetries, true, log)

        if (response != null) {
            response.results.each {
                String proteinMd5 = it.md5.toLowerCase()
                if (it.found) {
                    calculatedMatches[proteinMd5] = [:]
                    it.matches.each { matchObj ->
                        String library = matchObj.signature.signatureLibraryRelease.library
                        String appName = library.toLowerCase().replaceAll("[-\\s]", "")

                        if (applications.contains(appName)) {
                            matchObj = transformMatch(matchObj)
                            calculatedMatches[proteinMd5][matchObj.modelAccession] = matchObj
                        }
                    }
                } else {
                    def seq = sequences[proteinMd5]
                    noLookupFasta.append(">${proteinMd5}\n")
                    noLookupFasta.append("${seq}\n")
                }
            }
        } else {
            success = false
            break
        }
    }

    if (success) {
        def jsonMatches = JsonOutput.toJson(calculatedMatches)
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson(calculatedMatches))
        if (noLookupFasta.length() != 0) { new File(noLookupFastaPath.toString()).write(noLookupFasta.toString()) }
    } else {
        log.warn "An error occurred while querying the Matches API, analyses will be run locally"
        // when the connection fails, write out all sequences to "noLookup.fasta"
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        if (noLookupFasta.length() != 0) { new File(noLookupFastaPath.toString()).write(noLookupFasta.toString()) }
    }
}

String getMatchesApiUrl(matchesApiUrl, lookupServiceUrl, datadir, workflowManifest) {
    String _matchesApiUrl = matchesApiUrl ?: lookupServiceUrl
    // Get MLS metadata: api (version), release, release_date
    Map info = InterPro.httpRequest("${sanitizeURL(baseUrl)}/info", null, 0, true, log)
    if (info == null) {
        log.warn "An error occurred while querying the Matches API; analyses will be run locally"
        _matchesApiUrl = null
    } else {
        def apiVersion = info.api ?: "X.Y.Z"
        def majorVersion = apiVersion.split("\\.")[0]

        if (majorVersion != "0") {
            log.warn "${workflowManifest.name} ${workflowManifest.version}" +
                    " is not compatible with the Matches API at ${_matchesApiUrl};" +
                    " analyses will be run locally"
            _matchesApiUrl = null
        } else if (datadir) {
            def iprVersion = getDatabaseVersion("InterPro", datadir)
            if (iprVersion != info.release) {
                log.warn "InterPro version mismatch (local: ${iprVersion}, Matches API: ${info.release});" +
                        " pre-calulated matches will not be retrieved, and analyses will be run locally"
                _matchesApiUrl = null
            }
        }
    }
    return _matchesApiUrl
}

def Map transformMatch(Map match) {
    // * operator - spread contents of a map or collecion into another map or collection
    return [
        *            : match,
        "treegrafter": ["ancestralNodeID": match["annotationNode"]],
        "locations"  : match["locations"].collect { loc ->
            return [
                *          : loc,
                "hmmBounds": loc["hmmBounds"] ? Location.getReverseHmmBounds(loc["hmmBounds"]) : null,
                "fragments": loc["fragments"].collect { tranformFragment(it) },
                "sites"    : loc["sites"] ?: []
            ]
        },
    ]
}

def Map tranformFragment(Map fragment) {
    return [
        "start"   : fragment["start"],
        "end"     : fragment["end"],
        "dcStatus": fragment["type"]
    ]
}