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

    jsonSlurper.parse(seqMatches).each { sequence ->
        def seqMatches = []
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            memberDB = matchObj.signature.signatureLibraryRelease.library

            matchResult = [
                "signature": matchObj.signature,
                "locations": []
            ]
            if (matchObj.locations) {
                matchObj.locations.each { location ->
                    locationResult = [
                        "start": location.start,
                        "end": location.end,
                        "representative": location.representative
                    ]
                    hmmBounds = Output.convertHmmBound(location.hmmBounds)
                    switch (memberDB) {
                        case "cdd":
                            locationResult["evalue"] = matchObj.evalue
                            locationResult["score"] = matchObj.score
                            break
                        case "coils":
                            break
                        case "hamap":
                            locationResult["score"] = location.score
                            locationResult["alignment"] = location.targetAlignment ?: "Not available"
                            break
                        case "mobidblite":
                            locationResult["sequence-feature"] = location.sequenceFeature
                            break
                        case "panther":
                            locationResult["hmmStart"] = location.hmmStart
                            locationResult["hmmEnd"] = location.hmmEnd
                            locationResult["hmmLength"] = 0
                            locationResult["hmmBounds"] = hmmBounds
                            locationResult["envelopeStart"] = location.envelopeStart
                            locationResult["envelopeEnd"] = location.envelopeEnd
                            break
                       case "phobius":
                            locationResult["score"] = location.score
                            locationResult["prediction"] = location.prediction
                            locationResult["topology"] = location.topology
                            break
                        case "pirsf":
                            locationResult["evalue"] = location.evalue
                            locationResult["score"] = location.score
                            locationResult["hmmStart"] = location.start
                            locationResult["hmmEnd"] = location.end
                            locationResult["hmmLength"] = location.hmmLength
                            locationResult["hmmBounds"] = hmmBounds
                            locationResult["envelopeStart"] = location.envelopeStart
                            locationResult["envelopeEnd"] = location.envelopeEnd
                            break
                        case "prints":
                            locationResult["pvalue"] = location.pvalue
                            locationResult["score"] = location.score
                            locationResult["motifNumber"] = location.motifNumber
                            break
                        case "prositeprofiles":
                            locationResult["score"] = location.score
                            locationResult["alignment"] = location.targetAlignment
                            break
                        case "prositepatterns":
                            locationResult["cigarAlignment"] = location.cigarAlignment
                            locationResult["alignment"] = location.targetAlignment
                            locationResult["level"] = location.level
                            break
                        case ["pirsr", "sfld"]:
                            locationResult["evalue"] = location.evalue
                            locationResult["score"] = location.score
                            locationResult["hmmStart"] = location.hmmStart
                            locationResult["hmmEnd"] = location.hmmEnd
                            locationResult["hmmLength"] = location.hmmLength
                            locationResult["envelopeStart"] = location.envelopeStart
                            locationResult["envelopeEnd"] = location.envelopeEnd
                            break
                        case ["signalp", "signalp_euk"]:
                            locationResult["pvalue"] = location.pvalue
                            locationResult["cleavageStart"] = signalp.cleavageSiteStart
                            locationResult["cleavageEnd"] = signalp.cleavageSiteEnd
                            break
                        case "smart":
                            locationResult["evalue"] = location.evalue
                            locationResult["score"] = location.score
                            locationResult["hmmStart"] = location.hmmStart
                            locationResult["hmmEnd"] = location.hmmEnd
                            locationResult["hmmLength"] = location.hmmLength
                            locationResult["hmmBounds"] = hmmBounds
                            break
                        case "superfamily":
                            locationResult["hmmLength"] = location.hmmLength
                            matchResult['evalue'] = location.evalue
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
                    if (memberDB in ["pirsr", "sfld"]) {
                        locationResult["sites"] = location.sites ?: []
                    }
                    if (memberDB == "cdd") {
                        locationResult["sites"] = location.sites?.collect { site ->
                            site.remove("label")
                            site.remove("group")
                            site.remove("hmmStart")
                            site.remove("hmmEnd")
                            return site
                        } ?: []
                    }
                    locationResult["location-fragments"] = location.fragments
                    matchResult["locations"].add(locationResult)
                }
            }

            if (matchResult["locations"]) {
                if (!(memberDB in ["cdd", "coils", "hamap", "mobidblite", "phobius", "prositeprofiles", "prositepatterns", "prints", "signalp", "signalp_euk", "superfamily"])) {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["score"] = matchObj.score
                }
                matchResult["model-ac"] = matchObj.modelAccession.split("\\.")[0]
                if (memberDB == "sfld") {
                    matchResult["scope"] = null
                } else if (memberDB == "panther") {
                    matchResult["name"] = matchObj.treegrafter.subfamilyDescription
                    matchResult["accession"] = matchObj.treegrafter.subfamilyAccession
                    matchResult["model-ac"] = matchObj.treegrafter.subfamilyAccession
                    matchResult["goXRefs"] = matchObj.signature?.entry?.goXRefs ?: []
                    matchResult["proteinClass"] = matchObj.treegrafter.proteinClass
                    matchResult["graftPoint"] = matchObj.treegrafter.graftPoint
                } else if (memberDB == "prints") {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["graphscan"] = matchObj.graphScan
                } else if (memberDB in ["signalp", "signalp_euk"]) {
                    matchResult["orgType"] = matchObj.orgType
                }

                if (memberDB in ['cathfunfam', 'cathgene3d', 'panther']) {
                    name = matchObj.signature.description
                    description = matchObj.signature.name
                    matchObj.signature.name = name
                    matchObj.signature.description = description
                }
                seqMatches.add(matchResult)
            }
        }
        sequence["matches"] = seqMatches
        jsonOutput["results"].add(sequence)
    }

    def outputFilePath = "${outputPath}.ips6.json"
    def json = JsonOutput.toJson(jsonOutput)
    new File(outputFilePath.toString()).write(json)
}
