import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PAINT_ANNOTATIONS {
    label 'xref'

    // Retrieve PAINT annotations for Panther hits
    // calculated and pre-calc becuase they are not retrieved from the Match Lookup
    input:
    val paintAnnDir
    tuple val(meta), val(membersMatches)

    output:
    tuple val(meta), path("paint_annotation.json")

    exec:
    JsonSlurper jsonSlurper = new JsonSlurper()
    def matches = jsonSlurper.parse(membersMatches).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            Match matchObject = Match.fromMap(jsonMatch)
            if (matchObject.signature.signatureLibraryRelease.library == "panther") {
                String sigAcc = matchObject.signature.accession
                String paintAnnPath = "${paintAnnDir}/${sigAcc}.json"
                File paintAnnotationFile = new File(paintAnnPath.toString())
                if (paintAnnotationFile.exists()) {
                    def paintAnnotationsContent = jsonSlurper.parse(paintAnnotationFile)
                    String nodeId = matchObject.treegrafter.ancestralNodeID
                    def nodeData = paintAnnotationsContent[nodeId]
                    matchObject.treegrafter.proteinClass = nodeData[2]
                    matchObject.treegrafter.graftPoint = nodeData[3]
                }
            }
        }]
    }
    def outputFilePath = task.workDir.resolve("paint_annotation.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
