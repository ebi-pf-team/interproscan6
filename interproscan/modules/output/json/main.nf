import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process WRITE_JSON_OUTPUT {
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

    Map<String, List<String>> membersLocationFields = [
        "CDD": ["evalue-match", "score-match"],
        "COILS": [],
        "HAMAP": ["score", "targetAlignment"],
        "MobiDB Lite": ["sequence-feature"],
        "PANTHER": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "Phobius": [],
        "PIRSF": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "PIRSR": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "PRINTS": ["pvalue", "score", "motifNumber"],
        "PROSITE profiles": ["score", "targetAlignment"],
        "PROSITE patterns": ["cigarAlignment", "targetAlignment", "level"],
        "SFLD": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "SignalP_Euk": ["pvalue", "cleavageStart", "cleavageEnd"],
        "SignalP-Prok": ["pvalue", "cleavageStart", "cleavageEnd"],
        "SMART": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "SUPERFAMILY": ["hmmLength", "evalue"],
        "TMHMM": []
    ]
    List<String> otherMembersLocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]

    jsonSlurper.parse(seqMatches).each { sequence ->
        def seqMatches = []
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            String memberDB = matchObj.signature.signatureLibraryRelease.library ?: ""
            matchResult = [
                "signature": matchObj.signature,
                "locations": []
            ]
            if (matchObj.locations) {
                // Locations
                matchObj.locations.each { location ->
                    locationResult = [
                        "start": location.start,
                        "end": location.end,
                        "representative": location.representative,
                        "location-fragments": location.fragments
                    ]

                    hmmBounds = location.getHmmBounds(location.hmmBounds)
                    def fields = membersLocationFields.get(memberDB, otherMembersLocationFields)
                    fields.each { field ->
                        switch (field) {
                            case "targetAlignment":
                                locationResult["alignment"] = location.targetAlignment
                                break
                            case "cigarAlignment":
                                locationResult["cigarAlignment"] = location.cigarAlignment
                                break
                            case "cleavageStart":
                                locationResult["cleavageStart"] = matchObj.signalp.cleavageSiteStart
                                break
                            case "cleavageEnd":
                                locationResult["cleavageEnd"] = matchObj.signalp.cleavageSiteEnd
                                break
                            case "envelopeStart":
                                locationResult["envelopeStart"] = location.envelopeStart
                                break
                            case "envelopeEnd":
                                locationResult["envelopeEnd"] = location.envelopeEnd
                                break
                            case "evalue":
                                locationResult["evalue"] = location.evalue
                                break
                            case "evalue-match":
                                locationResult["evalue"] = matchObj.evalue
                                break
                            case "hmmStart":
                                locationResult["hmmStart"] = location.hmmStart
                                break
                            case "hmmEnd":
                                locationResult["hmmEnd"] = location.hmmEnd
                                break
                            case "hmmLength":
                                locationResult["hmmLength"] = location.hmmLength
                                break
                            case "hmmBounds":
                                locationResult["hmmBounds"] = hmmBounds
                                break
                            case "level":
                                locationResult["level"] = location.level
                                break
                            case "motifNumber":
                                locationResult["motifNumber"] = location.motifNumber
                                break
                            case "pvalue":
                                locationResult["pvalue"] = location.pvalue
                                break
                            case "score":
                                locationResult["score"] = location.score
                                break
                            case "score-match":
                                locationResult["score"] = matchObj.score
                                break
                            case "sequence-feature":
                                locationResult["sequence-feature"] = location.sequenceFeature
                                break
                            default:
                                println "Warning: Unknown fields '${field}'"
                                break
                        }
                    }

                    // Sites
                    if (memberDB in ["PIRSR", "SFLD"]) {
                        locationResult["sites"] = location.sites ?: []
                    }
                    if (memberDB == "CDD") {
                        locationResult["sites"] = location.sites?.collect { site ->
                            site.remove("label")
                            site.remove("group")
                            site.remove("hmmStart")
                            site.remove("hmmEnd")
                            return site
                        } ?: []
                    }
                    matchResult["locations"].add(locationResult)
                }
            }

            if (matchResult["locations"]) {
                // Match level fields
                if (memberDB in ["AntiFam", "FunFam", "CATH-Gene3D", "NCBIfam", "PANTHER", "Pfam", "PIRSF", "PIRSR", "SFLD", "SMART", "TMHMM"]) {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["score"] = matchObj.score
                }
                matchResult["model-ac"] = matchObj.modelAccession

                switch (memberDB) {
                    case "SFLD":
                        matchResult["scope"] = null
                        break
                    case "PANTHER":
                        matchResult["name"] = matchObj.treegrafter.subfamilyDescription
                        matchResult["accession"] = matchObj.treegrafter.subfamilyAccession
                        matchResult["model-ac"] = matchObj.treegrafter.subfamilyAccession
                        matchResult["goXRefs"] = matchObj.signature?.entry?.goXRefs ?: []
                        matchResult["proteinClass"] = matchObj.treegrafter.proteinClass
                        matchResult["graftPoint"] = matchObj.treegrafter.graftPoint
                        break
                    case "PRINTS":
                        matchResult["evalue"] = matchObj.evalue
                        matchResult["graphscan"] = matchObj.graphScan
                        break
                    case ["SignalP_Euk", "SignalP-Prok"]:
                        matchResult["orgType"] = matchObj.signalp.orgType
                }
                seqMatches.add(matchResult)
            }
        }
        sequence["matches"] = seqMatches
        jsonOutput["results"].add(sequence)
    }

    def outputFilePath = "${outputPath}.ips6.json"
    def json = JsonOutput.prettyPrint(JsonOutput.toJson(jsonOutput))
    new File(outputFilePath.toString()).write(json)
}
