import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.util.regex.Pattern

HMM_BOUND_PATTERN = [
    "[]": "COMPLETE",
    "[.": "N_TERMINAL_COMPLETE",
    ".]": "C_TERMINAL_COMPLETE",
    "..": "INCOMPLETE"
]

process JSON_OUTPUT {
    label 'write_output'

    input:
    val seqMatches
    val outputPath
    val ips6Version

    exec:
    def jsonSlurper = new JsonSlurper()
    def jsonOutput = [:]
    jsonOutput["interproscan-version"] = ips6Version
    jsonOutput["results"] = []

    jsonSlurper.parse(seqMatches).each { sequence ->
        def seqMatches = []
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            memberDB = matchObj.signature.signatureLibraryRelease.library
            matchResult = [
                "signature": matchObj.signature,
                "locations": []
            ]
            if (memberDB in ["PHOBIUS", "SUPERFAMILY"]) {
                match_data['locations'].each { location ->
                    def locationResult = [
                        "start"          : location.start,
                        "end"            : location.end,
                        "representative" : location.representative,
                        "location-fragments": location.locationFragments
                    ]

                    if (memberDB == "PHOBIUS") {
                        matchResult["model-ac"] = matchObj.modelAccession
                    } else if (memberDB == "SUPERFAMILY") {
                        locationResult['evalue'] = location.evalue
                        locationResult["hmmLength"] = location.hmmLength
                        matchResult["locations"] = location_info
                        matchResult["evalue"] = matchObj.evalue
                        matchResult["model-ac"] = matchObj.modelAccession
                    } else {
                        locationResult["sequence-feature"] = matchObj.sequenceFeature
                    }

                    matchResult["locations"].add(locationResult)
                }
            } else {
                if (matchObj.locations) {
                    matchObj.locations.each { location ->
                        locationResult = [
                            "start": location.start,
                            "end": location.end,
                            "representative": location.representative
                        ]
                        hmmBounds = HMM_BOUND_PATTERN[location.hmmBounds]
                        switch (memberDB) {
                            case "CDD":
                                locationResult["evalue"] = location.evalue
                                locationResult["score"] = location.score
                                break
                            case "HAMAP":
                                locationResult["score"] = location.score
                                locationResult["alignment"] = location.alignment
                                break
                            case "MOBIDB_LITE":
                                locationResult["sequence-feature"] = location.sequenceFeature
                                break
                            case "PANTHER":
                                locationResult["hmmStart"] = location.hmmStart
                                locationResult["hmmEnd"] = location.hmmEnd
                                locationResult["hmmLength"] = 0
                                locationResult["hmmBounds"] = hmmBounds
                                locationResult["envelopeStart"] = location.envelopeStart
                                locationResult["envelopeEnd"] = location.envelopeEnd
                                break
                            case "PIRSF":
                                locationResult["evalue"] = location.evalue
                                locationResult["score"] = location.score
                                locationResult["hmmStart"] = location.start
                                locationResult["hmmEnd"] = location.end
                                locationResult["hmmLength"] = location.hmmLength
                                locationResult["hmmBounds"] = hmmBounds
                                locationResult["envelopeStart"] = location.envelopeStart
                                locationResult["envelopeEnd"] = location.envelopeEnd
                                break
                            case "PRINTS":
                                locationResult["pvalue"] = location.pvalue
                                locationResult["score"] = location.score
                                locationResult["motifNumber"] = location.motifNumber
                                break
                            case "PROSITE_PROFILES":
                                locationResult["score"] = location.score
                                locationResult["alignment"] = location.alignment
                                break
                            case "PROSITE_PATTERNS":
                                locationResult["cigarAlignment"] = location.cigarAlignment
                                locationResult["alignment"] = location.alignment
                                locationResult["level"] = location.level
                                break
                            case ["PIRSR", "SFLD"]:
                                locationResult["evalue"] = location.evalue
                                locationResult["score"] = location.score
                                locationResult["hmmStart"] = location.hmmStart
                                locationResult["hmmEnd"] = location.hmmEnd
                                locationResult["hmmLength"] = location.hmmLength
                                locationResult["envelopeStart"] = location.envelopeStart
                                locationResult["envelopeEnd"] = location.envelopeEnd
                                break
                            case ["SIGNALP", "SIGNALP_EUK"]:
                                locationResult["pvalue"] = location.pvalue
                                locationResult["cleavageStart"] = location.cleavageStart
                                locationResult["cleavageEnd"] = location.cleavageEnd
                                break
                            case "SMART":
                                locationResult["evalue"] = location.evalue
                                locationResult["score"] = location.score
                                locationResult["hmmStart"] = location.hmmStart
                                locationResult["hmmEnd"] = location.hmmEnd
                                locationResult["hmmLength"] = location.hmmLength
                                locationResult["hmmBounds"] = hmmBounds
                                break
                            default:
                                locationResult["evalue"] = location.evalue
                                locationResult["score"] = location.score
                                locationResult["hmmStart"] = location.hmmStart
                                locationResult["hmmEnd"] = location.hmmEnd
                                locationResult["hmmLength"] = location.hmmLength
                                locationResult["hmmBounds"] = hmmBounds
                                locationResult["envelopeStart"] = location.envelopeStart
                                locationResult["envelopeEnd"] = location.envelopeEnd
                        }
                        if (memberDB in ["CDD", "PIRSR", "SFLD"]) {
                            locationResult["sites"] = location.sites ?: []
                        }
                        locationResult["location-fragments"] = location.fragments

                        matchResult["locations"].add(locationResult)
                    }
                }

                if (!(memberDB in ["CDD", "COILS", "HAMAP", "MOBIDB_LITE", "PHOBIUS", "PIRSR", "PROSITE_PROFILES", "PROSITE_PATTERNS", "PRINTS", "SIGNALP", "SIGNALP_EUK"])) {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["score"] = matchObj.score
                }
                matchResult["model-ac"] = matchObj.modelAccession
                if (memberDB == "SFLD") {
                    matchResult["scope"] = null
                } else if (memberDB == "PANTHER") {
                    matchResult["name"] = matchObj.entry.subfamily_description
                    matchResult["accession"] =matchObj.modelAccession
                    matchResult["goXRefs"] = matchObj.goXrefs
                    signature["description"] = null
                    signature["name"] = matchObj.entry.description
                    matchResult["proteinClass"] = matchObj.proteinClass
                    matchResult["graftPoint"] = matchObj.graftPoint
                } else if (memberDB == "PRINTS") {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["graphscan"] = matchObj.graphscan
                } else if (memberDB in ["SIGNALP", "SIGNALP_EUK"]) {
                    matchResult["orgType"] = matchObj.orgType
                }
            }
            seqMatches.add(matchResult)
        }
        sequence["matches"] = seqMatches
        jsonOutput["results"].add(sequence)
    }

    def outputFilePath = "${outputPath}.ips6.json"
    def json = JsonOutput.toJson(jsonOutput)
    new File(outputFilePath.toString()).write(json)
}
