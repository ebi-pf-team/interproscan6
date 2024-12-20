import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'small'

    input:
    tuple val(index), val(fasta), val(json)
    val applications
    val chunkSize
    val host
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
                jsonResponse.each { md5, matches ->
                    seqId = sequences[md5].id
                    if (matches != null) {
                        calculatedMatches[seqId] = [:]
                        matches.each { match ->
                            Match matchObj = Match.fromMap(match)
                            memberDB = matchObj.signature.signatureLibraryRelease.library
                            stdMemberDB = memberDB.toLowerCase().replaceAll("[-\\s]", "")
                            if (applications.contains(stdMemberDB)) {
                                modelAccession = matchObj.signature.accession
                                matchObj.modelAccession = modelAccession
                                calculatedMatches[seqId][modelAccession] = matchObj
                            }
                        }
                    } else {
                        def seq = sequences[md5]
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
                log.error "Unable to connect to the match lookup service. Max retries reached. Running analysis locally"
                exceededRetries = true
                break
            }
            log.warn "Could not connect to the match lookup service. Retrying connection."
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
        new File(calculatedMatchesPath.toString()).write(JsonOutput.toJson([:]))
        new File(fasta.toString()).copyTo(new File(noLookupFastaPath.toString()))
        jsonFile.copyTo(new File(noLookupMapPath.toString()))
    }
}
