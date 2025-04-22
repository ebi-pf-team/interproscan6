import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process PREPARE_LOOKUP {
    executor 'local'

    input:
    val _url
    val db_releases
    val interproscan_version  // major.minor iprscan version number
    val workflow_manifest

    output:
    val matchesApiUrl

    exec:
    String _matchesApiUrl = _url  // reassign to avoid 'variable' already used error when logging
    // Get MLS metadata: api (version), release, release_date
    Map info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(_matchesApiUrl)}/info".toString(), null, 0, true)
    if (info == null) {
        log.warn "An error occurred while querying the Matches API; analyses will be run locally"
        matchesApiUrl = null
    } else {
        def apiVersion = info.api ?: "X.Y.Z"
        def majorVersion = apiVersion.split("\\.")[0]
        if (majorVersion != "0") {
            log.warn "${workflow_manifest.name} ${workflow_manifest.version}" +
                    " is not compatible with the Matches API at ${_matchesApiUrl};" +
                    " analyses will be run locally"
            matchesApiUrl = null
        } else if (db_releases) {  // can be null if we don't need data for the selected apps (e.g. mobidblite)
            if (db_releases["interpro"]["version"] != info.release) {
                log.warn "The local InterPro version does not match the match API release (Local: ${db_releases['interpro']}, Matches API: ${info.release}).\n" +
                        "Pre-calculated matches will not be retrieved, and analyses will be run locally"
                matchesApiUrl = null
            }
        }
    }
    matchesApiUrl = _matchesApiUrl
    return matchesApiUrl
}

process LOOKUP_MATCHES {
    executor 'local'

    input:
    tuple val(index), val(fasta), val(applications), val(url), val(chunkSize), val(maxRetries)

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

    String baseUrl = HTTPRequest.sanitizeURL(url.toString())
    boolean success = true
    for (chunk in chunks) {
        String data = JsonOutput.toJson([md5: chunk])
        def response = HTTPRequest.fetch("${baseUrl}/matches", data, maxRetries, true)

        if (response != null) {
            response.results.each {
                String proteinMd5 = it.md5.toUpperCase() // ensure it matches the local seq Db case
                if (it.found) {
                    calculatedMatches[proteinMd5] = [:]
                    it.matches.each { matchMap ->
                        String library = matchMap.signature.signatureLibraryRelease.library
                        if (library == "MobiDB Lite") {
                            matchMap.signature.signatureLibraryRelease.library = "MobiDB-lite"
                        }
                        
                        String appName = library.toLowerCase().replaceAll("[-\\s]", "")

                        if (applications.contains(appName)) {
                            matchMap = transformMatch(matchMap)
                            calculatedMatches[proteinMd5][matchMap.modelAccession] = matchMap
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

def Map transformMatch(Map match) {
    // * operator - spread contents of a map or collecion into another map or collection
    return [
        *            : match,
        "treegrafter": ["ancestralNodeID": match["annotationNode"]],
        "locations"  : match["locations"].collect { loc ->
            return [
                *          : loc,
                "hmmBounds": loc["hmmBounds"] ? getReverseHmmBounds(loc["hmmBounds"]) : null,
                "fragments": loc["fragments"].collect { tranformFragment(it) },
                "sites"    : loc["sites"] ?: []
            ]
        },
    ]
}

def getReverseHmmBounds(hmmBounds) {
    return [
        "COMPLETE"            : "[]",
        "N_TERMINAL_COMPLETE" : "[.",
        "C_TERMINAL_COMPLETE" : ".]",
        "INCOMPLETE"          : ".."
    ][hmmBounds]
}

def Map tranformFragment(Map fragment) {
    return [
        "start"   : fragment["start"],
        "end"     : fragment["end"],
        "dcStatus": fragment["type"]
    ]
}
