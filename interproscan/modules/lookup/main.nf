import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'local'

    input:
    tuple val(index), val(fasta), val(json)
    val applications
    val url
    val chunkSize
    val maxRetries

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), path("noLookup.json"), optional: true

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")
    def noLookupMapPath = task.workDir.resolve("noLookup.json")

    def calculatedMatches = [:]
    def noLookupFasta = new StringBuilder()
    def noLookupMap = [:]

    def jsonSlurper = new JsonSlurper()
    def jsonFile = new File(json.toString())
    def sequences = jsonSlurper.parse(jsonFile)
        .collectEntries{ seqId, obj ->
            if (obj instanceof List) { // nucleotide sequences case
                obj.collectEntries { seq ->
                    [(seq.md5): FastaSequence.fromMap(seq)]
                }
            } else {
                [(obj.md5): FastaSequence.fromMap(obj)]
            }
        }

    def md5List = sequences.keySet().toList()
    def chunks = md5List.collate(chunkSize)

    String baseUrl = sanitizeURL(url.toString())
    boolean success = true
    for (chunk in chunks) {
        String data = JsonOutput.toJson([md5: chunk])
        def response = httpRequest("${baseUrl}/matches", data, maxRetries, true)

        if (response != null) {
            response.results.each {
                seqId = sequences[it.md5].id
                if (it.found) {
                    calculatedMatches[seqId] = [:]
                    it.matches.each { matchObj ->
                        String library = matchObj.signature.signatureLibraryRelease.library
                        String appName = library.toLowerCase().replaceAll("[-\\s]", "")

                        if (applications.contains(appName)) {
                            matchObj = transformMatch(matchObj)
                            calculatedMatches[seqId][matchObj.modelAccession] = matchObj
                        }
                    }
                } else {
                    def seq = sequences[it.md5]
                    noLookupMap[seqId] = seq
                    noLookupFasta.append(">${seqId} ${seq.description}\n")
                    noLookupFasta.append("${seq.sequence}\n")
                }
            }
        } else {
            success = false
            break
        }
    }

    if (success) {
        def jsonMatches = JsonOutput.toJson(calculatedMatches)
        new File(calculatedMatchesPath.toString()).write(jsonMatches)

        if (!noLookupMap.isEmpty()) {
            def jsonSequences = JsonOutput.toJson(noLookupMap)
            new File(noLookupMapPath.toString()).write(jsonSequences)
            new File(noLookupFastaPath.toString()).write(noLookupFasta.toString())
        }
    } else {
        log.warn "An error occurred while querying the Matches API, analyses will be run locally"
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        new File(fasta.toString()).copyTo(new File(noLookupFastaPath.toString()))
        jsonFile.copyTo(new File(noLookupMapPath.toString()))
    }
}

def String sanitizeURL(String url) {
    return url.replaceAll(/\/+$/, '')
}

def Map getInfo(String baseUrl) {
    return httpRequest("${sanitizeURL(baseUrl)}/info", null, 0, true)
}

def httpRequest(String urlString, String data, int maxRetries, boolean verbose) {
    int attempts = 0
    boolean isPost = data != null && data.length() > 0
    HttpURLConnection connection = null
    while (attempts < maxRetries + 1) {
        attempts++

        try {
            URL url = new URL(urlString)
            connection = (HttpURLConnection) url.openConnection()
            connection.connectTimeout = 5000
            connection.readTimeout = 5000
            if (isPost) {
                connection.doOutput = true
                connection.requestMethod = "POST"
                connection.setRequestProperty("Content-Type", "application/json")
                connection.setRequestProperty("Accept", "application/json")
                connection.outputStream.withWriter("UTF-8") { writer ->
                    writer << data
                }
                connection.outputStream.flush()
            } else {
                connection.requestMethod = "GET"
            }

            int responseCode = connection.responseCode
            if (responseCode >= 200 && responseCode < 300) {
                String responseText = connection.inputStream.getText("UTF-8")
                return new JsonSlurper().parseText(responseText)
            } else {
                if (verbose) {
                    def errorMsg = connection.errorStream ? connection.errorStream.getText("UTF-8") : ""
                    log.warn "Received HTTP ${responseCode} for ${urlString}"
                }
                return null
            }
        } catch (java.net.ConnectException | java.net.UnknownHostException e) {
            if (verbose) {
                log.warn "Connection error for ${urlString}: ${e.message}"
            }
            return null
        } catch (java.net.SocketTimeoutException e) {
            if (verbose) {
                log.warn "Timeout error for ${urlString}; attempt ${attempts} / ${maxRetries + 1}"
            }
        } catch (Exception e) {
            if (verbose) {
                log.warn "Unexpected error for ${urlString}: ${e.message}"
            }
            return null
        } finally {
            connection?.disconnect()
        }
    }
}

def Map transformMatch(Map match) {
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