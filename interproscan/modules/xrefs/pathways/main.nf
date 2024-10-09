import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PATHWAYS {
    label 'xref'

    input:
    val ipr2paJson
    val paInfoJson
    tuple val(meta), val(membersMatches)

    output:
    tuple val(meta), path("matches2pa.json")

    exec:
    def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
    ]

    JsonSlurper jsonSlurper = new JsonSlurper()
    def ipr2pa = jsonSlurper.parseText(ipr2paJson.text)
    def paInfo = jsonSlurper.parseText(paInfoJson.text)
    def matches = jsonSlurper.parse(membersMatches).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            def matchObject = Match.fromMap(jsonMatch)
            println "matchObject: ${matchObject.modelAccession}"
            if (matchObject.signature.entry) {
                def interproKey = matchObject.signature.entry.accession
                println "interproKey: ${interproKey}"
                try {
                    def paIds = ipr2pa[interproKey]
                    def paTerms = paIds.collect { paId ->
                        [
                            "name": paInfo[paId][1],
                            "databaseName": PA_PATTERN[paInfo[paId][0]],
                            "id": paId
                        ]
                    }
                    matchObject.signature.entry.paTerms = paTerms
                } catch (Exception e) {
                    // pass
                }
            }
        }]
    }
    def outputFilePath = task.workDir.resolve("matches2pa.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
