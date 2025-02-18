import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'local'

    input:
    tuple val(index), val(fasta)
    val applications
    val chunkSize
    val host
    val maxRetries

    output:
    tuple val(index), path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), path("noLookup.json"), optional: true

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def calculatedMatches = [:]

    def noLookupFastaPath = task.workDir.resolve("noLookup.fasta")
    def noLookupFasta = new StringBuilder()

    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]
    def md5List = sequences.keySet().toList()
    def chunks = md5List.collate(chunkSize)

    int attempt = 0
    boolean success = false
    boolean exceededRetries = false
    while (attempt < maxRetries && !success) {
        try {
            chunks.each { chunk ->
                def requestBody = JsonOutput.toJson([md5: chunk])
                def connection = new URL(host).openConnection()
                connection.requestMethod = 'POST'
                connection.doOutput = true
                connection.setRequestProperty('Content-Type', 'application/json')
                connection.outputStream.withWriter('UTF-8') { writer ->
                    writer << requestBody
                }
                def response = connection.inputStream.text

                def jsonResponse = new JsonSlurper().parseText(response)
                jsonResponse.results.each {
                    seqId = sequences[it.md5].id
                    if (it.found) {
                        calculatedMatches[seqId] = [:]
                        it.matches.each { match ->
                            Match matchObj = Match.fromMap(match)
                            memberDB = match.signature.signatureLibraryRelease.library.toLowerCase().replaceAll("[-\\s]", "")
                            if (applications.contains(memberDB)) {
                                matchObj.modelAccession = matchObj.signature.accession
                                calculatedMatches[seqId][modelAccession] = matchObj
                            }
                        }
                    } else {
                        def seq = sequences[it.md5]
                        noLookupMap[seqId] = seq
                        noLookupFasta.append(">${seqId} ${seq.description}\n")
                        noLookupFasta.append("${seq.sequence}\n")
                    }
                }
            }
            success = true
        } catch (Exception e) {
            attempt++
            if (attempt >= maxRetries) {
                log.error "Unable to retrieve pre-calculated matches. Running analyses locally."
                exceededRetries = true
                break
            }
            log.warn "An error occurred when retrieving pre-calculated matches. Retrying."
        }
    }

    if (success) {
        def jsonMatches = JsonOutput.toJson(calculatedMatches)
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson(calculatedMatches))
    } else {
        // when the connection fails, write out all sequences to "noLookup.fasta"
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        new File(fasta.toString()).copyTo(new File(noLookupFastaPath.toString()))
    }
}
