import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'mls'

    input:
    tuple val(index), val(fasta), val(json)
    val appl
    val chunkSize
    val host

    output:
    path("calculatedMatches.json")
    path("noLookupSeq.json")

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def noLookupSeqPath = task.workDir.resolve("noLookupSeq.json")

    def jsonSlurper = new JsonSlurper()
    def jsonFile = new File(json.toString())
    sequences = jsonSlurper.parse(jsonFile)

    def md5List = []
    sequences.each { md5, seq_info ->
        md5List << md5
    }
    def matchesResult = [:]
    def noLookupSeq = [:]
    def chunks = md5List.collate(chunkSize)
    chunks.each { chunk ->
        def url = "${host}"
        def requestBody = JsonOutput.toJson([md5: chunk])
        def connection = new URL(url).openConnection()
        connection.requestMethod = 'POST'
        connection.doOutput = true
        connection.setRequestProperty('Content-Type', 'application/json')
        connection.outputStream.withWriter('UTF-8') { writer ->
            writer << requestBody
        }
        def response = connection.inputStream.text
        def jsonResponse = new JsonSlurper().parseText(response)
        jsonResponse.each { entry ->
            if (entry.value != null) {
                matchesResult[entry.key] = entry.value
            } else {
                noLookupSeq = sequences.findAll { md5, seq_info ->
                    md5 == entry.key
                }
            }
        }
    }
    def jsonMatches = JsonOutput.toJson(matchesResult)
    def jsonSequences = JsonOutput.toJson(noLookupSeq)
    new File(calculatedMatchesPath.toString()).write(jsonMatches)
    new File(noLookupSeqPath.toString()).write(jsonSequences)
}
