def check_matches_api(applications, api_url, _version) {
    def appls_in_api     = []
    def appls_not_in_api = []
    def sanitized_url = uk.ac.ebi.interpro.HTTPRequest.sanitizeURL(api_url)
    def url = "${sanitized_url}/info".toString()
    def info = uk.ac.ebi.interpro.HTTPRequest.fetch(url, null, 0)
    if (info) {
        def api_version = info.api ?: "X.Y.Z"
        def major_version = api_version.split("\\.")[0]
        if (major_version != "0") {
            log.warn "This version of InterProScan is not compatible with the Matches API at ${sanitized_url}; pre-calculated annotations will not be retrieved."
        } else {
            if (info.analyses) {
                def all_appls_in_api = info.analyses*.name.collect { n -> 
                    n.toLowerCase().replaceAll("[-\\s]", "") 
                }
                applications.each { appl -> 
                    if (all_appls_in_api.contains(appl)) {
                        appls_in_api << appl
                    } else {
                        appls_not_in_api << appl
                    }
                }
            } else {
                log.warn "Could not retrieve the list of analyses available in the Matches API; pre-calculated annotations will not be retrieved."
            }
        }
    } else {
        log.warn "An error occurred while querying the Matches API; pre-calculated annotations will not be retrieved."
    }

    if (appls_in_api.isEmpty()) {
        appls_not_in_api = applications.clone() as List<String>
    }

    return [appls_in_api, appls_not_in_api]
}


process GET_MATCHES {
    maxForks 1
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fasta)
    val applications
    val api_url
    val chunk_size
    val max_retries

    output:
    tuple val(meta), path("matches.json"),                  emit: json
    tuple val(meta), path("unknown.fasta"), optional: true, emit: fasta

    exec:
    def matches = [:]
    def fasta_sb = new StringBuilder()
    def sequences = uk.ac.ebi.interpro.FastaFile.parse(fasta)  // [md5: sequence]
    def md5s = sequences.keySet().toList().sort()
    def request_chunk_size = chunk_size > 100 ? 100 : chunk_size // API set max to 100
    def chunks = md5s.collate(request_chunk_size)
    def base_url = uk.ac.ebi.interpro.HTTPRequest.sanitizeURL(api_url.toString())
    def response = null
    def success = chunks.every { chunk ->
        def data = groovy.json.JsonOutput.toJson([md5: chunk])
        response = uk.ac.ebi.interpro.HTTPRequest.fetch("${base_url}/matches", data, max_retries, true)
        if (response == null) {
            return false
        }

        response.results.each { item ->
            def seq_md5 = item.md5.toUpperCase()
            if (item.found) {
                matches[seq_md5] = [:]
                item.matches.each { match ->
                    def library = match.signature.signatureLibraryRelease.library
                    def app_name = library.toLowerCase().replaceAll("[-\\s]", "")
                    if (applications.contains(app_name)) {
                        match = transformMatch(match, sequences[seq_md5])
                        matches[seq_md5][match.modelAccession] = match
                    }
                }
            } else {
                def seq = sequences[seq_md5]
                fasta_sb.append(">${seq_md5}\n")
                fasta_sb.append("${seq}\n")
            }
        }
        return true
    }

    def json_file = task.workDir.resolve("matches.json")
    def fasta_file = task.workDir.resolve("unknown.fasta")

    if (success) {
        def jf = new com.fasterxml.jackson.core.JsonFactory()
        json_file.withWriter { writer ->
            def gen = jf.createGenerator(writer)
            new com.fasterxml.jackson.databind.ObjectMapper().writeValue(gen, matches)
            gen.close()
        }
        if (fasta_sb.length() != 0) {
            fasta_file.text = fasta_sb.toString()
        }
    } else {
        log.warn "An error occurred while querying the Matches API -- '${response}'"
        json_file.text = groovy.json.JsonOutput.toJson([:])
        fasta_file.text = fasta.text
    }
}

def transformMatch(match, seq) {
    def transformedMatch = [:] + match
    transformedMatch["modelAccession"] = match["model-ac"]
    transformedMatch["treegrafter"] = ["ancestralNodeID": match["ancestralNode"]]
    transformedMatch["locations"] = match["locations"].collect { loc ->
        def transformedLocation = [:] + loc
        transformedLocation["sequenceFeature"] = loc["sequence-feature"]
        transformedLocation["hmmBounds"] = loc["hmmBounds"] ? getReverseHmmBounds(loc["hmmBounds"]) : null
        transformedLocation["fragments"] = loc["location-fragments"].collect { f -> tranformFragment(f) }
        transformedLocation["sites"] = loc["sites"] ?: []
        transformedLocation["targetAlignment"] = loc["cigarAlignment"] ? decodeAlignment(loc["cigarAlignment"], seq, loc["start"]) : null
        return transformedLocation
    }
    return transformedMatch
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
        if (op == 'M' || op == '=' || op == 'X') {
            targetAlign << sequence.substring(index, index + len)
            index += len
        } else if (op == 'I') {
            targetAlign << sequence.substring(index, index + len).toLowerCase()
            index += len
        } else if (op == 'D') {
            targetAlign << '-' * len
        } else {
            throw new IllegalArgumentException("Unsupported CIGAR operation: $op")
        }
    }
    return targetAlign.toString()
}

def tranformFragment(fragment) {
    return [
        "start"   : fragment["start"],
        "end"     : fragment["end"],
        "dcStatus": fragment["dc-status"]
    ]
}
