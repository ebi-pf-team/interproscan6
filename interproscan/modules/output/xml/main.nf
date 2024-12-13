import groovy.json.JsonSlurper
import groovy.xml.MarkupBuilder
import java.time.format.DateTimeFormatter
import java.time.LocalDate
import java.io.StringWriter
import java.util.regex.Pattern

process WRITE_XML_OUTPUT {
    label 'write_output'

    input:
    val matches
    val outputPath
    val nucleic
    val ips6Version

    exec:
    def NT_SEQ_ID_PATTERN = Pattern.compile(/^orf\d+\s+source=(.*)\s+coords=(\d+)\.\.(\d+)\s+.+frame=(\d+)\s+desc=(.*)$/)
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    // set the correct encoding so symbols are formatted correctly in the final output
    xml.setEscapeAttributes(false) // Prevent escaping attributes
    xml.setEscapeText(false)       // Prevent escaping text

    def jsonSlurper = new JsonSlurper()
    def jsonData = jsonSlurper.parse(matches)
    def processedNT = []

    String matchType = nucleic ? "nucleotide-sequence-matches" : "protein-matches"
    xml."$matchType"("interproscan-version": ips6Version){
        jsonData.each { seqData ->
            if (nucleic) {
                def ntSequenceMD5 = seqData.translatedFrom[0]["md5"]
                "nucleotide-sequence" {
                    if (!processedNT.contains(ntSequenceMD5)) {
                        processedNT << ntSequenceMD5
                        sequence(md5: ntSequenceMD5, seqData.translatedFrom[0]["sequence"])
                        seqData.translatedFrom.each { crossRef ->
                            xref(id: crossRef.id, name: "${crossRef.id} ${crossRef.description}")
                        }
                    }
                    def ntMatch = NT_SEQ_ID_PATTERN.matcher(seqData.xref[0].name)
                    assert ntMatch.matches()
                    start   : ntMatch.group(2) as int,
                    end     : ntMatch.group(3) as int,
                    strand  : (ntMatch.group(4) as int) < 4 ? "SENSE" : "ANTISENSE",
                    orf (end: end, start: start, strand: strand) {
                        protein {
                            sequence(md5: seqData["md5"], seqData["sequence"])
                            matches {
                                processMatches(seqData["matches"], xml)
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
            } else {
                protein {
                    sequence(md5: seqData["md5"], seqData["sequence"])
                    matches {
                        processMatches(seqData["matches"], xml)
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
    }

    def outputFilePath = "${outputPath}.ips6.xml"
    new File(outputFilePath).text = writer.toString()
}

def processMatches(matches, xml) {
    List<String> hmmer3Members = ["AntiFam", "CATH-Gene3D", "FunFam", "HAMAP", "NCBIfam", "Pfam", "PIRSF", "PIRSR", "SFLD", "SUPERFAMILY"]
    List<String> hmmer3LocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]
    Map<String, List<String>> memberLocationFields = [
        "CDD": ["match-evalue", "match-score"],
        "COILS": [],
        "HAMAP": ["score"],
        "MobiDB-lite": ["sequence-feature"],
        "PANTHER": ["hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"],
        "Phobius": [],
        "PRINTS": ["motifNumber", "pvalue", "score"],
        "PROSITE patterns": ["level"],
        "PROSITE profiles": ["score"],
        "SignalP": ["score"],
        "SMART": ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds"],
        "SUPERFAMILY": ["evalue", "hmmLength"]
    ]

    matches.each { modelAcc, matchMap ->
        Match matchObj = Match.fromMap(matchMap)
        def matchNodeName
        def memberDb = matchObj.signature.signatureLibraryRelease.library
        if (hmmer3Members.contains(memberDb)) {
            matchNodeName = "hmmer3"
        } else if (memberDb == "SMART") {
            matchNodeName = "hmmer2"
        } else {
            matchNodeName = memberDb
        }

        def matchAttributes = [:]
        if (hmmer3Members.findAll { it != "SUPERFAMILY"}.contains(memberDb) || memberDb == "SMART") {
            matchAttributes.evalue = matchObj.evalue
            matchAttributes.score = matchObj.score
        } else if (memberDb == "PANTHER") {
            matchAttributes.ac = matchObj.treegrafter.subfamilyAccession
            matchAttributes.evalue = matchObj.evalue
            matchAttributes."graft-point" = matchObj.treegrafter.graftPoint
            matchAttributes.name = matchObj.signature.name
            matchAttributes.score = matchObj.score
        } else if (memberDb == "PRINTS") {
            matchAttributes.evalue = matchObj.evalue
            matchAttributes.graftscan = matchObj.graphScan
        }

        xml."$matchNodeName-match"(matchAttributes) {
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

            "model-ac"(memberDb == "PANTHER" ? matchObj.treegrafter.subfamilyAccession : matchObj.modelAccession)

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
                            if (memberDb in ["HAMAP", "PROSITE patterns", "PROSITE profiles"]) {
                                alignment(loc.targetAlignment ?: "")
                            }
                            if (loc.sites) {
                                "sites" {
                                    loc.sites.each { siteObj ->
                                        "$matchNodeName-site"(description: siteObj.description, numLocations: siteObj.numLocations) {
                                            if(siteObj.group){ group(siteObj.group) }
                                            if(siteObj.label){ label(siteObj.label) }
                                            if (memberDb != "CDD") {
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
