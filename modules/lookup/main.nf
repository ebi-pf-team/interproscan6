import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process PREPARE_LOOKUP {
    /* A Simple process to check API and InterPro version compatibility
    Retain as a process so that this process and the LOOKUP subworkflow wait for the
    channels to be ready before determining if the API is available */
    label    'tiny'
    executor 'local'

    input:
    val matches_api_apps
    val api_interpro_version
    val db_releases
    val url

    output:
    val api_url

    exec:
    _url = url // reassign to avoid variable already declared error
    if (db_releases["interpro"]["version"] != api_interpro_version) {
            log.warn "The local InterPro version (${db_releases['interpro']}) does not match the Matches API release (${api_interpro_version}). Pre-calculated matches will not be retrieved and analyses will run locally."
            _url = null
    }
    api_url = _url
}

process LOOKUP_MATCHES {
    maxForks 1
    label    'tiny'
    executor 'local'

    input:
    tuple val(index), val(fasta), val(applications), val(url), val(chunkSize), val(maxRetries)

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), optional: true

    exec:
    def seqFasta = fasta.toString()  // reassign to avoid variable already declared error

    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")

    def calculatedMatches = [:]
    def noLookupFasta = new StringBuilder()
    Map<String, String> sequences = FastaFile.parse(seqFasta)  // [md5: sequence]
    def md5List = sequences.keySet().toList().sort()
    def chunks = md5List.collate(chunkSize)

    String baseUrl = HTTPRequest.sanitizeURL(url.toString())
    boolean success = true

    for (chunk in chunks) {
        data = JsonOutput.toJson([md5: chunk])
        response = HTTPRequest.fetch("${baseUrl}/matches", data, maxRetries, true)

        if (response != null) {
            response.results.each {
                String proteinMd5 = it.md5.toUpperCase()
                if (it.found) {
                    calculatedMatches[proteinMd5] = [:]
                    it.matches.each { matchMap ->
                        String library = matchMap.signature.signatureLibraryRelease.library
                        if (library == "MobiDB Lite") {
                            matchMap.signature.signatureLibraryRelease.library = "MobiDB-lite"
                        }
                        String appName = library.toLowerCase().replaceAll("[-\\s]", "")
                        if (applications.contains(appName)) {
                            matchMap = transformMatch(matchMap, sequences[proteinMd5])
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
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson(calculatedMatches))
        if (noLookupFasta.length() != 0) {
            new File(noLookupFastaPath.toString()).write(noLookupFasta.toString())
        }
    } else {
        log.warn "An error occurred while querying the Matches API, analyses will be run locally"
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        new File(noLookupFastaPath.toString()).write(new File(fasta.toString()).text)
    }
}

def Map transformMatch(Map match, String seq) {
    // * operator - spread contents of a map or collecion into another map or collection
    return [
        *            : match,
        "treegrafter": ["ancestralNodeID": match["annotationNode"]],
        "locations"  : match["locations"].collect { loc ->
            return [
                *          : loc,
                "hmmBounds": loc["hmmBounds"] ? getReverseHmmBounds(loc["hmmBounds"]) : null,
                "fragments": loc["fragments"].collect { tranformFragment(it) },
                "sites"    : loc["sites"] ?: [],
                "targetAlignment": loc["cigarAlignment"] ? decodeAlignment(loc["cigarAlignment"], seq, loc["start"]) : null
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

def decodeAlignment(cigarAlignment, sequence, startIndex) {
    def targetAlign = new StringBuilder()
    def index = startIndex - 1 // convert from 1-based numbering to 0-based numbering
    def matcher = (cigarAlignment =~ /(\d+)([MID=X])/)
    matcher.each { match ->
        def len = match[1].toInteger()
        def op = match[2]
        switch(op) {
            case 'M':
            case '=':
            case 'X':
                targetAlign << sequence.substring(index, index + len)
                index += len
                break
            case 'I':
                targetAlign << sequence.substring(index, index + len).toLowerCase()
                index += len
                break
            case 'D':
                targetAlign << '-' * len
                break
            default:
                throw new IllegalArgumentException("Unsupported CIGAR operation: $op")
        }
    }
    return targetAlign.toString()
}

def Map tranformFragment(Map fragment) {
    return [
        "start"   : fragment["start"],
        "end"     : fragment["end"],
        "dcStatus": fragment["type"]
    ]
}
