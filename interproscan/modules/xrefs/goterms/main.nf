import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process GOTERMS {
    label 'xref'

    input:
    val ipr2goJson
    val goInfoJson
    tuple val(meta), val(membersMatches)

    output:
    tuple val(meta), path("matches2go.json")

    exec:
    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    JsonSlurper jsonSlurper = new JsonSlurper()
    def ipr2go = jsonSlurper.parse(ipr2goJson)
    def goInfo = jsonSlurper.parse(goInfoJson)
    def matches = jsonSlurper.parse(membersMatches).collectEntries { seqId, jsonMatches ->
        [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
            Match matchObject = Match.fromMap(jsonMatch)
            if (matchObject.signature.entry) {
                String interproKey = matchObject.signature.entry.accession
                try {
                    def goIds = ipr2go[interproKey]
                    def goTerms = goIds.collect { goId ->
                        [
                            "name": goInfo[goId][1],
                            "databaseName": GO_PATTERN[goInfo[goId][0]],
                            "id": goId
                        ]
                    }
                    matchObject.signature.entry.goTerms = goTerms
                } catch (Exception e) {
                    // pass
                }
            }
        }]
    }
    def outputFilePath = task.workDir.resolve("matches2go.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
