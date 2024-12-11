import groovy.json.JsonSlurper
import java.net.URL
import groovy.json.JsonOutput

process LOOKUP_MATCHES {
    label 'small'

    input:
    tuple val(index), val(fasta), val(json)
    val appl
    val chunkSize
    val host

    output:
    path("calculatedMatches.json")
    tuple val(index), path("noLookup.fasta"), path("noLookup.json")

    exec:
    def calculatedMatchesPath = task.workDir.resolve("calculatedMatches.json")
    def noLookupFasta = task.workDir.resolve("noLookup.fasta")
    def noLookupSeqPath = task.workDir.resolve("noLookup.json")

    def jsonSlurper = new JsonSlurper()
    def jsonFile = new File(json.toString())
    def sequences = jsonSlurper.parse(jsonFile)
        .collectEntries{ seqId, obj ->
            if (obj instanceof List) { // nucleotide sequences case
                obj.collectEntries { orf ->
                    [(orf.md5): FastaSequence.fromMap(orf)]
                }
            } else {
                [(obj.md5): FastaSequence.fromMap(obj)]
            }
        }

    def md5List = sequences.keySet().toList()
//     md5List << "5FE1059FDE57D6E61C5343CAB0C502C8"
    def matchesResult = [:]
    def noLookupSeq = [:]
    noLookupMD5 = []
    def chunks = md5List.collate(chunkSize)
    chunks.each { chunk ->
        def requestBody = JsonOutput.toJson([md5: chunk])
        def connection = new URL("${host}").openConnection()
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
                noLookupMD5 << entry.key
            }
        }
    }

    noLookupSeq = sequences.findAll { seqId, seq_info ->
        noLookupMD5.contains(seq_info.md5)
    }
    new File(noLookupFasta.toString()).withWriter { writer ->
        noLookupSeq.each { seqId, seq ->
            writer.writeLine(">${seq.id} ${seq.description}")
            writer.writeLine(seq.sequence)
        }
    }
    def jsonMatches = JsonOutput.toJson(matchesResult)
    def jsonSequences = JsonOutput.toJson(noLookupSeq)
    new File(calculatedMatchesPath.toString()).write(jsonMatches)
    new File(noLookupSeqPath.toString()).write(jsonSequences)
}
