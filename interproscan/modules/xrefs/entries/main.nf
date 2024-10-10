import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import com.fasterxml.jackson.databind.ObjectMapper

process ENTRIES {
    label 'xref'

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
            "version": entries['databases'][memberDB]
        ]
        def matches = jsonSlurper.parse(matchesPath).collectEntries { seqId, jsonMatches ->
            [(seqId): jsonMatches.collectEntries { matchId, jsonMatch ->
                Match matchObject = Match.fromMap(jsonMatch)
                def entriesInfo = entries['entries']
                String accId = matchObject.modelAccession
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

                if (memberDB == "panther") {
                    String accSubfamily = matchObject.treegrafter.subfamilyAccession
                    matchObject.treegrafter.subfamilyAccession = accSubfamily
                    if (entriesInfo[accSubfamily]) {
                        matchObject.treegrafter.subfamilyName = entrieInfo[accSubfamily]["name"]
//                         Just waiting finish output step to have sure we don't use this infos
//                         matchObject.treegrafter.subfamilyDesc = entriesInfo[accSubfamily]["description"]
//                         matchObject.treegrafter.subfamilyType = entriesInfo[accSubfamily]["type"]
                    }
                }

                [(matchId): matchObject]
            }]
        }
        def outputFilePath = task.workDir.resolve("matches2entries.json")
        def json = JsonOutput.toJson(matches)
        new File(outputFilePath.toString()).write(json)
    }
}
