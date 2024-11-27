import groovy.json.JsonSlurper
import groovy.json.JsonOutput

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

    Map<String, List<String>> membersLocationFields = [
        "cdd": [],
        "coils": [],
        "hamap": ["score", "targetAlignment"],
        "mobidblite": ["sequence-feature"],
        "panther": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "phobius": [],
        "pirsf": ["evalue", "score", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "pirsr": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "prints": ["pvalue", "score", "motifNumber"],
        "prositeprofiles": ["score", "targetAlignment"],
        "prositepatterns": ["cigarAlignment", "targetAlignment", "level"],
        "sfld": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "signalp": ["pvalue", "cleavageStart", "cleavageEnd"],
        "signalp_euk": ["pvalue", "cleavageStart", "cleavageEnd"],
        "smart": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "superfamily": ["hmmLength"]
    ]
    List<String> otherMembersLocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]

    def boundsMapping = [
        "[]"  : "COMPLETE",
        "[."  : "N_TERMINAL_COMPLETE",
        ".]"  : "C_TERMINAL_COMPLETE",
        ".."  : "INCOMPLETE"
    ]

    jsonSlurper.parse(seqMatches).each { sequence ->
        def seqMatches = []
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            String rawMemberDB = matchObj.signature.signatureLibraryRelease.library
            String memberDB = rawMemberDB.toLowerCase().replace("-", "").replace(" ", "")
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

                    // field value exceptions
                    if (memberDB == "pirsf") {
                        locationResult["hmmStart"] = location.start
                        locationResult["hmmEnd"] = location.end
                    } else if (memberDB == "cdd") {
                        locationResult["evalue"] = matchObj.evalue
                        locationResult["score"] = matchObj.score
                    } else if (memberDB == "superfamily") {
                        matchResult["evalue"] = location.evalue
                    }

                    hmmBounds = boundsMapping[location.hmmBounds]

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
                            case "sequence-feature":
                                locationResult["sequence-feature"] = location.sequenceFeature
                                break
                            default:
                                println "Warning: Unknown fields '${field}'"
                                break
                        }
                    }

                    // Sites
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
                    matchResult["locations"].add(locationResult)
                }
            }

            if (matchResult["locations"]) {
                // Match level fields
                if (memberDB in ["antifam", "cathfunfam", "cathgene3d", "ncbifam", "panther", "pfam", "pirsf", "pirsr", "sfld", "smart", "tmhmm"]) {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["score"] = matchObj.score
                }
                matchResult["model-ac"] = memberDB == "cathfunfam" ? matchObj.modelAccession : matchObj.modelAccession.split("\\.")[0]

                switch (memberDB) {
                    case "sfld":
                        matchResult["scope"] = null
                        break
                    case "panther":
                        matchResult["name"] = matchObj.treegrafter.subfamilyDescription
                        matchResult["accession"] = matchObj.treegrafter.subfamilyAccession
                        matchResult["model-ac"] = matchObj.treegrafter.subfamilyAccession
                        matchResult["goXRefs"] = matchObj.signature?.entry?.goXRefs ?: []
                        matchResult["proteinClass"] = matchObj.treegrafter.proteinClass
                        matchResult["graftPoint"] = matchObj.treegrafter.graftPoint
                        break
                    case "prints":
                        matchResult["evalue"] = matchObj.evalue
                        matchResult["graphscan"] = matchObj.graphScan
                        break
                    case ["signalp", "signalp_euk"]:
                        matchResult["orgType"] = matchObj.signalp.orgType
                }

                if (memberDB in ['cathfunfam', 'cathgene3d', 'panther', 'superfamily']) {
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
    def json = JsonOutput.prettyPrint(JsonOutput.toJson(jsonOutput))
    new File(outputFilePath.toString()).write(json)
}
