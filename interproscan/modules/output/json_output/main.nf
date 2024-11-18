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

    def locationCases = [
        "cdd": [],
        "coils": [],
        "hamap": ["score", "alignment"],
        "mobidblite": ["sequence-feature"],
        "panther": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "phobius": [],
        "pirsf": ["evalue", "score", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "pirsr": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "prints": ["pvalue", "score", "motifNumber"],
        "prositeprofiles": ["score", "alignment"],
        "prositepatterns": ["cigarAlignment", "alignment", "level"],
        "sfld": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "envelopeStart", "envelopeEnd"],
        "signalp": ["pvalue", "cleavageStart", "cleavageEnd"],
        "signalp_euk": ["pvalue", "cleavageStart", "cleavageEnd"],
        "smart": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "superfamily": ["hmmLength", "evalue"],
        "others": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]
    ]

    jsonSlurper.parse(seqMatches).each { sequence ->
        def seqMatches = []
        sequence["matches"].each { matchId, match ->
            Match matchObj = Match.fromMap(match)
            String memberDB = InterProScan.standardiseMemberDB(matchObj.signature.signatureLibraryRelease.library)
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

                    hmmBounds = switch (location.hmmBounds) {
                        case "[]":
                            return "COMPLETE"
                        case "[.":
                            return "N_TERMINAL_COMPLETE"
                        case ".]":
                            return "C_TERMINAL_COMPLETE"
                        case "..":
                            return "INCOMPLETE"
                    }

                    // field value exceptions
                    if (memberDB = "pirsf") {
                        locationResult["hmmStart"] = location.start
                        locationResult["hmmEnd"] = location.end
                    } else if (memberDB = "cdd") {
                        locationResult["evalue"] = matchObj.evalue
                        locationResult["score"] = matchObj.score
                    }

                    specialCases.each { db, fields ->
                        if (db == memberDB || (db instanceof List && memberDB in db)) {
                            fields.each { field ->
                                switch (field) {
                                    case "alignment":
                                        locationResult["alignment"] = location.alignment
                                    case "cigarAlignment":
                                        locationResult["cigarAlignment"] = location.cigarAlignment
                                    case "cleavageStart":
                                        locationResult["cleavageStart"] = matchObj.signalp.cleavageSiteStart
                                    case "cleavageEnd":
                                        locationResult["cleavageEnd"] = matchObj.signalp.cleavageSiteEnd
                                    case "envelopeStart":
                                        locationResult["envelopeStart"] = location.envelopeStart
                                    case "envelopeEnd":
                                        locationResult["envelopeEnd"] = location.envelopeEnd
                                    case "evalue":
                                        locationResult["evalue"] = location.evalue
                                    case "hmmStart":
                                        locationResult["hmmStart"] = location.hmmStart
                                    case "hmmEnd":
                                        locationResult["hmmEnd"] = location.hmmEnd
                                    case "hmmLength":
                                        locationResult["hmmLength"] = location.hmmLength
                                    case "hmmBounds":
                                        locationResult["hmmBounds"] = hmmBounds
                                    case "level":
                                        locationResult["level"] = location.level
                                    case "motifNumber":
                                        locationResult["motifNumber"] = location.motifNumber
                                    case "pvalue":
                                        locationResult["pvalue"] = location.pvalue
                                    case "score":
                                        locationResult["score"] = location.score
                                    case "sequence-feature":
                                        locationResult["sequence-feature"] = location.sequenceFeature

                                }
                            }
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
                if (memberDB in ["antifam", "funfam", "gene3d", "ncbifam", "panther", "pfam", "pirsf", "pirsr", "sfld", "smart", "tmhmm"]) {
                    matchResult["evalue"] = matchObj.evalue
                    matchResult["score"] = matchObj.score
                }
                matchResult["model-ac"] = memberDB == "FunFam" ? matchObj.modelAccession : matchObj.modelAccession.split("\\.")[0]

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
                    case in ["signalp", "signalp_euk"]:
                        matchResult["orgType"] = matchObj.signalp.orgType
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
