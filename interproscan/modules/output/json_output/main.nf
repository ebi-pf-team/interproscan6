import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import java.util.regex.Pattern

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

    def matches = jsonSlurper.parse(seqMatches).each { sequence ->
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            memberDB = matchObj.signature.signatureLibraryRelease.library
            finalMatch = [
                "signature": matchObj.signature
            ]
            if (memberDB in ["PHOBIUS", "SUPERFAMILY"]) {
                match_data['locations'].each { location ->
                    def locationInfo = [
                        "start"          : matchObj.locations.start,
                        "end"            : matchObj.locations.end,
                        "representative" : matchObj.locations.representative,
                        "location-fragments": matchObj.locations.locationFragments
                    ]

                    if (memberDB == "PHOBIUS") {
                        finalMatch["locations"] = locationInfo
                        finalMatch["model-ac"] = matchObj.modelAccession
                    } else if (memberDB == "SUPERFAMILY") {
                        locationInfo['evalue'] = matchObj.locations.evalue
                        locationInfo["hmmLength"] = matchObj.locations.hmmLength
                        finalMatch["locations"] = location_info
                        finalMatch["evalue"] = matchObj.evalue
                        finalMatch["model-ac"] = matchObj.modelAccession
                    } else {
                        locationInfo["sequence-feature"] = matchObj.sequenceFeature
                        finalMatch["locations"] = location_info
                    }
                }
            } else {
                if (matchObj.locations) {
                    matchObj.locations.each { location ->
                        info = [
                            "start": location.start,
                            "end": location.end,
                            "representative": location.representative
                        ]
                        switch (memberDB) {
                            case "CDD":
                                info["evalue"] = location.evalue
                                info["score"] = location.score
                                break
                            case "HAMAP":
                                info["score"] = location.score
                                info["alignment"] = location.alignment
                                break
                            case "MOBIDB_LITE":
                                info["sequence-feature"] = location.sequenceFeature
                                break
                            case "PANTHER":
                                info["hmmStart"] = location.hmmStart
                                info["hmmEnd"] = location.hmmEnd
                                info["hmmLength"] = 0
                                info["hmmBounds"] = location.hmmBounds
                                info["envelopeStart"] = location.envelopeStart
                                info["envelopeEnd"] = location.envelopeEnd
                                break
                            case "PIRSF":
                                info["evalue"] = location.evalue
                                info["score"] = location.score
                                info["hmmStart"] = location.start
                                info["hmmEnd"] = location.end
                                info["hmmLength"] = location.hmmLength
                                info["hmmBounds"] = location.hmmBounds
                                info["envelopeStart"] = location.envelopeStart
                                info["envelopeEnd"] = location.envelopeEnd
                                break
                            case "PRINTS":
                                info["pvalue"] = location.pvalue
                                info["score"] = location.score
                                info["motifNumber"] = location.motifNumber
                                break
                            case "PROSITE_PROFILES":
                                info["score"] = location.score
                                info["alignment"] = location.alignment
                                break
                            case "PROSITE_PATTERNS":
                                info["cigarAlignment"] = location.cigarAlignment
                                info["alignment"] = location.alignment
                                info["level"] = location.level
                                break
                            case ["PIRSR", "SFLD"]:
                                info["evalue"] = location.evalue
                                info["score"] = location.score
                                info["hmmStart"] = location.hmmStart
                                info["hmmEnd"] = location.hmmEnd
                                info["hmmLength"] = location.hmmLength
                                info["envelopeStart"] = location.envelopeStart
                                info["envelopeEnd"] = location.envelopeEnd
                                break
                            case ["SIGNALP", "SIGNALP_EUK"]:
                                info["pvalue"] = location.pvalue
                                info["cleavageStart"] = location.cleavageStart
                                info["cleavageEnd"] = location.cleavageEnd
                                break
                            case "SMART":
                                info["evalue"] = location.evalue
                                info["score"] = location.score
                                info["hmmStart"] = location.hmmStart
                                info["hmmEnd"] = location.hmmEnd
                                info["hmmLength"] = location.hmmLength
                                info["hmmBounds"] = location.hmmBounds
                                break
                            default:
                                info["evalue"] = location.evalue
                                info["score"] = location.score
                                info["hmmStart"] = location.hmmStart
                                info["hmmEnd"] = location.hmmEnd
                                info["hmmLength"] = location.hmmLength
                                info["hmmBounds"] = location.hmmBounds
                                info["envelopeStart"] = location.envelopeStart
                                info["envelopeEnd"] = location.envelopeEnd
                        }

                        if (memberDB in ["CDD", "PIRSR", "SFLD"]) {
                            info["sites"] = location.sites ?: []
                        }

                        info["location-fragments"] = location.fragments

                        finalMatch["locations"] = info
                    }

                    if (!(memberDB in ["CDD", "COILS", "HAMAP", "MOBIDB_LITE", "PHOBIUS", "PIRSR", "PROSITE_PROFILES", "PROSITE_PATTERNS", "PRINTS", "SIGNALP", "SIGNALP_EUK"])) {
                        finalMatch["evalue"] = matchObj.evalue
                        finalMatch["score"] = matchObj.score
                    }

                    match["model-ac"] = matchObj.modelAccession
                    if (memberDB == "SFLD") {
                        match["scope"] = null
                    } else if (memberDB == "PANTHER") {
                        match["name"] = matchObj.entry.subfamily_description
                        match["accession"] =matchObj.modelAccession
                        match["goXRefs"] = matchObj.goXrefs
                        signature["description"] = null
                        signature["name"] = matchObj.entry.description
                        match["proteinClass"] = matchObj.proteinClass
                        match["graftPoint"] = matchObj.graftPoint
                    } else if (memberDB == "PRINTS") {
                        match["evalue"] = matchObj.evalue
                        match["graphscan"] = matchObj.graphscan
                    } else if (memberDB in ["SIGNALP", "SIGNALP_EUK"]) {
                        match["orgType"] = matchObj.orgType
                    }

                    finalMatch["matches"] = match
                }
            }
        }
        jsonOutput["results"].add(sequence)
    }

    def outputFilePath = "${outputPath}.ips6.json"
    def json = JsonOutput.toJson(jsonOutput)
    new File(outputFilePath.toString()).write(json)
}
