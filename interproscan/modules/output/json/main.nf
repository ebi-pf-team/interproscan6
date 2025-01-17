import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.JsonParser
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.core.JsonEncoding
import com.fasterxml.jackson.core.JsonToken
import com.fasterxml.jackson.databind.SerializationFeature
import java.util.regex.Pattern

process WRITE_JSON_OUTPUT {
    label 'write_output'

    input:
    val sequenceMatches
    val outputPath
    val nucleic
    val ips6Version

    exec:
    ObjectMapper mapper = new ObjectMapper()
    mapper.enable(SerializationFeature.INDENT_OUTPUT)

    def outputFilePath = "${outputPath}.ips6.json"
    JsonGenerator generator = mapper.getFactory().createGenerator(new File(outputFilePath.toString()), JsonEncoding.UTF8)
    generator.writeStartObject()
    generator.writeStringField("interproscan-version", ips6Version)
    generator.writeFieldName("results")
    generator.writeStartArray()

    JsonParser parser = mapper.getFactory().createParser(new File(sequenceMatches.toString()))
    parser.nextToken()

    def NT_SEQ_ID_PATTERN = Pattern.compile(/^orf\d+\s+source=(.*)\s+coords=(\d+)\.\.(\d+)\s+.+frame=(\d+)\s+desc=(.*)$/)
    def nucleicResults = [:]

    Map<String, List<String>> membersLocationFields = [
        "CDD": ["evalue-match", "score-match"],
        "COILS": [],
        "HAMAP": ["score", "targetAlignment"],
        "MobiDB-lite": ["sequence-feature"],
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
        "DeepTMHMM": []
    ]
    List<String> otherMembersLocationFields = ["evalue", "score", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart", "envelopeEnd"]

    while (parser.nextToken() != JsonToken.END_ARRAY) {
        def sequence = parser.readValueAs(Map)
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

        if (nucleic) {
            def nucleicSeqMd5 = sequence.translatedFrom[0].md5  // nucleic sequence md5 - same for all ORFs
            if (!nucleicResults.containsKey(nucleicSeqMd5)) {
                nucleicResults[nucleicSeqMd5] = [
                    sequence          : sequence.translatedFrom[0].sequence,
                    md5               : nucleicSeqMd5,
                    crossReferences   : [],
                    openReadingFrames : []
                ]
                sequence.translatedFrom.each { seq ->
                    nucleicResults[nucleicSeqMd5].crossReferences << [
                        name: "${seq.id} ${seq.description}",
                        id  : seq.id
                    ]
                }
            }

            def ntMatch = NT_SEQ_ID_PATTERN.matcher(sequence.xref[0].name)
            assert ntMatch.matches()
            nucleicResults[nucleicSeqMd5].openReadingFrames << [
                start   : ntMatch.group(2) as int,
                end     : ntMatch.group(3) as int,
                strand  : (ntMatch.group(4) as int) < 4 ? "SENSE" : "ANTISENSE",
                protein : [
                    sequence : sequence.sequence,
                    md5      : sequence.md5,
                    matches  : sequence.matches,
                    xref     : sequence.xref
                ]
            ]
        }

        if (nucleic) {
            nucleicResults.values().each { result ->
                generator.writeObject(result)
            }
        } else {
            generator.writeObject(sequence)
        }
    }

    generator.writeEndArray()
    generator.writeEndObject()
    parser.close()
    generator.close()
}
