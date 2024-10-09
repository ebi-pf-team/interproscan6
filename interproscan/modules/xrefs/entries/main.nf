import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import com.fasterxml.jackson.databind.ObjectMapper

process ENTRIES {
    input:
    tuple val(meta), val(membersMatches)
    val entriesJson

    output:
    tuple val(meta), path("matches2entries.json")

    exec:
    def objectMapper = new ObjectMapper()
    def entries = objectMapper.readValue(entriesJson, Map)

    JsonSlurper jsonSlurper = new JsonSlurper()
    matchesEntries = membersMatches.each { matchesPath  ->
        memberDB = matchesPath.toString().split("/").last().split("\\.")[0]
        def sigLibRelease = [
            "library": memberDB,
            "version": ""
        ]
        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                def matchObject = Match.fromMap(jsonMatch)
                def entriesInfo = entries['entries']
//                 def acc_id = match_key.split("\\.")[0]
                def accId = matchObject.modelAccession
                def entry = entriesInfo[accId] ?: entriesInfo[matchKey]
                entryData = null
                if (entry) {
                    def interproKey = entry['integrated']
                    def entryInfo = entriesInfo.get(interproKey)
                    if (entryInfo) {
                        entryData = [
                            "accession": interproKey,
                            "name": entryInfo["name"],
                            "description": entryInfo["description"],
                            "type": entryInfo["type"]
                        ]
                    }
                    def representativeFlag = null
                    if (entry['representative']) {
                        representative = new Tuple(entry['representative']["type"], entry['representative']["rank"])
                    }
                    matchObject.representativeFlag = representative
                }

                def sigData = [
                    "accession": entry["accession"],
                    "name": entry["name"],
                    "description": entry["description"],
                    "signatureLibraryRelease": sigLibRelease,
                    "entry": entryData
                ]
                Signature signatureObject = Signature.fromMap(sigData)
                matchObject.signature = signatureObject

                // TODO: add panther subfamily info on match object
                if (memberDB == "panther") {
                    def accSubfamily = data["accession"]
                    subfamilyName = entrieInfo[accSubfamily]["name"] ?: ""
                    subfamilyDesc = entriesInfo[accSubfamily]["description"] ?: ""
                    subfamilyType = entriesInfo[accSubfamily]["type"] ?: ""
                }

                [(matchId): matchObject]
            }]
        }
        def outputFilePath = task.workDir.resolve("matches2entries.json")
        def json = JsonOutput.toJson(matches)
        new File(outputFilePath.toString()).write(json)
    }
}
