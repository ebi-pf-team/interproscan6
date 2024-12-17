import groovy.json.JsonSlurper
import groovy.xml.MarkupBuilder
import java.time.format.DateTimeFormatter
import java.time.LocalDate
import java.io.StringWriter

process WRITE_XML_OUTPUT {
    label 'small'

    input:
    val matches
    val outputPath
    val ips6Version

    exec:
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    // set the correct encoding so symbols are formatted correctly in the final output
    xml.setEscapeAttributes(false) // Prevent escaping attributes
    xml.setEscapeText(false)       // Prevent escaping text

    def jsonSlurper = new JsonSlurper()
    def jsonData = jsonSlurper.parse(matches)

    List<String> hmmer3Members = ["antifam", "cathgene3d", "cathfunfam", "hamap", "ncbifam", "pfam", "pirsf", "pirsr", "sfld", "superfamily"]
    List<String> hmmer3LocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]
    Map<String, List<String>> memberLocationFields = [
        "cdd": ["match-evalue", "match-score"],
        "coils": [],
        "hamap": ["score"],
        "mobidblite": ["sequence-feature"],
        "panther": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "phobius": [],
        "prints": ["motifNumber", "pvalue", "score"],
        "prositepatterns": ["level"],
        "prositeprofiles": ["score"],
        "signalp": ["score"],
        "signalp-euk": ["score"],
        "smart": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "superfamily": ["evalue", "hmmLength"],
        "DeepTMHMM": []
    ]

    xml."protein-matches"("interproscan-version": ips6Version){
        jsonData.each { seqData ->
            protein {
                sequence(md5: seqData["md5"], seqData["sequence"])
                if (seqData["matches"]) {
                    matches {
                        seqData["matches"].each { modelAcc, matchMap ->
                            Match matchObj = Match.fromMap(matchMap)
                            def matchNodeName
                            def memberDb = matchObj.signature.signatureLibraryRelease.library
                            if (hmmer3Members.contains(memberDb)) {
                                matchNodeName = "hmmer3"
                            } else if (memberDb == "smart") {
                                matchNodeName = "hmmer2"
                            } else {
                                matchNodeName = memberDb
                            }

                            def matchAttributes = [:]
                            if (hmmer3Members.findAll { it != "superfamily"}.contains(memberDb) || memberDb == "smart") {
                                matchAttributes.evalue = matchObj.evalue
                                matchAttributes.score = matchObj.score
                            } else if (memberDb == "panther") {
                                matchAttributes.ac = matchObj.treegrafter.subfamilyAccession
                                matchAttributes.evalue = matchObj.evalue
                                matchAttributes."graft-point" = matchObj.treegrafter.graftPoint
                                matchAttributes.name = matchObj.signature.name
                                matchAttributes.score = matchObj.score
                            } else if (memberDb == "prints") {
                                matchAttributes.evalue = matchObj.evalue
                                matchAttributes.graftscan = matchObj.graphScan
                            }
                            "$matchNodeName-match"(matchAttributes) {
                                def signatureAttributes = [ac: matchObj.signature.accession]
                                if (matchObj.signature.name){
                                    signatureAttributes.name = matchObj.signature.name
                                }
                                if (matchObj.signature.description) {
                                    signatureAttributes.desc = matchObj.signature.description
                                }
                                signature(signatureAttributes) {
                                    signatureLibraryRelease {
                                        library(matchObj.signature.signatureLibraryRelease.library)
                                        version(matchObj.signature.signatureLibraryRelease.version)
                                    }
                                    if (matchObj.signature.entry) {
                                        matchObj.signature.entry.each { entryObj ->
                                            entry(
                                                ac: entryObj.accession,
                                                desc: entryObj.description ?: "-",
                                                name: entryObj.name ?: "-",
                                                type: entryObj.type ?: "-"
                                            ) {
                                                if (entryObj.goXRefs) {
                                                    entryObj.goXRefs.each { goXrefObj ->
                                                        "go-xref"(
                                                            category: goXrefObj.category,
                                                            db: goXrefObj.databaseName,
                                                            id: goXrefObj.id,
                                                            name: goXrefObj.name
                                                        )
                                                    }
                                                }
                                                if (entryObj.pathwayXRefs) {
                                                    entryObj.pathwayXRefs.each { pathwayObj ->
                                                        "pathway-xref"(
                                                            db: pathwayObj.databaseName,
                                                            id: pathwayObj.id,
                                                            name: pathwayObj.name
                                                        )
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                "model-ac"(memberDb == "panther" ? matchObj.treegrafter.subfamilyAccession : matchObj.modelAccession)

                                if (matchObj.locations) {
                                    def fields = memberLocationFields.get(memberDb, hmmer3LocationFields)
                                    locations {
                                        matchObj.locations.each { loc ->
                                            def locationAttributes = getLocationAttributes(loc, fields, matchObj)
                                            location(locationAttributes) {
                                                if (loc.fragments) {
                                                    "$matchNodeName-fragment" {
                                                        loc.fragments.each { frag ->
                                                            start(frag.start)
                                                            end(frag.end)
                                                            dcStatus(frag.dcStatus)
                                                        }
                                                    }
                                                }
                                                if (memberDb in ["hamap", "prositepatterns", "prositeprofiles"]) {
                                                    alignment(loc.targetAlignment ?: "")
                                                }
                                                if (loc.sites) {
                                                    "sites" {
                                                        loc.sites.each { siteObj ->
                                                            "$matchNodeName-site"(description: siteObj.description, numLocations: siteObj.numLocations) {
                                                                if(siteObj.group){ group(siteObj.group) }
                                                                if(siteObj.label){ label(siteObj.label) }
                                                                if (memberDb != "cdd") {
                                                                    hmmStart(siteObj.hmmStart)
                                                                    hmmEnd(siteObj.hmmEnd)
                                                                }
                                                                "site-locations" {
                                                                    siteObj.siteLocations.each { siteLoc ->
                                                                        "site-location"(
                                                                            residue: siteLoc.residue,
                                                                            start: siteLoc.start,
                                                                            end: siteLoc.end
                                                                        )
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                xrefs {
                    seqData.xref.each { ref ->
                        xref {
                            name(ref["name"])
                            id(ref["id"])
                        }
                    }
                }
            }
        }
    }

    def outputFilePath = "${outputPath}.ips6.xml"
    new File(outputFilePath).text = writer.toString()
}

def getLocationAttributes(Location location, List<String> memberFields, Match matchObj) {
    def locationAttributes = [
        start: location.start,
        end: location.end,
        representative: location.representative
    ]
    memberFields.each { field ->
        switch (field) {
            case "evalue":
                locationAttributes.evalue = location.evalue
                break
            case "score":
                locationAttributes.score = location.score
                break
            case "match-evalue":
                locationAttributes.evalue = matchObj.evalue
                break
            case "match-score":
                locationAttributes.score = matchObj.score
                break
            case "motifNumber":
                locationAttributes.motifNumber = location.motifNumber
                break
            case "pvalue":
                locationAttributes.pvalue = location.pvalue
                break
            case "level":
                locationAttributes.level = location.level
                break
            case "sequence-feature":
                locationAttributes["sequence-feature"] = location.sequenceFeature
                break
            case "hmmStart":
                locationAttributes["hmm-start"] = location.hmmStart
                break
            case "hmmEnd":
                locationAttributes["hmm-end"] = location.hmmEnd
                break
            case "hmmLength":
                locationAttributes["hmm-length"] = location.hmmLength
                break
            case "hmmBounds":
                locationAttributes["hmm-bounds"] = location.getHmmBounds(location.hmmBounds)
                break
            case "envelopeStart":
                locationAttributes["env-start"] = location.envelopeStart
                break
            case "envelopeEnd":
                locationAttributes["env-end"] = location.envelopeEnd
                break
            default:
                println "Warning: Unknown fields '${field}'"
                break
        }
    }
    return locationAttributes
}
