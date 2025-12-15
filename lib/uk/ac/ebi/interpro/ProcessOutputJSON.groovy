package uk.ac.ebi.interpro
import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.core.io.SerializedString
import com.fasterxml.jackson.core.util.DefaultPrettyPrinter
import com.fasterxml.jackson.core.util.DefaultIndenter
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import groovy.json.JsonException
import java.util.regex.Pattern
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.SeqDB


class ProcessOutputJSON {
    static void run(List<String> inputPaths, String databasePath, Map dbReleases, boolean isNucleic, String iprscanVersion, boolean useJsonLines, String outputPath) {
        ObjectMapper mapper = new ObjectMapper()
        SeqDB db = new SeqDB(databasePath)
        
        streamJson(outputPath) { JsonGenerator generator ->
            if (useJsonLines) {
                generator.setRootValueSeparator(new SerializedString(''))
                Set<String> seenNucleicMd5s = new HashSet<>()
                inputPaths.each { inputPath ->
                    Map proteins = mapper.readValue(new File(inputPath), Map)

                    if (isNucleic) {
                        def (nucleicToProteinMd5, ntSeqDataMap, orfDataMap) = 
                            db.retrieveAllNucleicSequenceData(proteins.keySet() as List)

                        nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                            if (!seenNucleicMd5s.contains(nucleicMd5)) {
                                generator.writeStartObject()
                                generator.writeStringField("interproscan-version", iprscanVersion)
                                generator.writeStringField("interpro-version", dbReleases?.interpro?.version)
                                generator.writeFieldName("results")
                                generator.writeStartArray()
                                def seqData = ntSeqDataMap[nucleicMd5]
                                def protMd5ToOrfs = orfDataMap[nucleicMd5]
                                writeNucleic(nucleicMd5, proteinMd5s, proteins, generator, seqData, protMd5ToOrfs)
                                generator.writeEndArray()
                                generator.writeEndObject()
                                generator.writeRaw('\n')
                                seenNucleicMd5s.add(nucleicMd5)
                            }
                        }
                    } else {
                        def allMd5s = proteins.keySet() as List
                        def md5ToSeqData = [:]
                        allMd5s.collate(1000).each { batch ->
                            def result = db.proteinMd5ToProteinSeqs(batch)
                            md5ToSeqData.putAll(result)
                        }
                        
                        proteins.each { String proteinMd5, Map proteinMatches ->
                            def seqData = md5ToSeqData[proteinMd5]
                            generator.writeStartObject()
                            generator.writeStringField("interproscan-version", iprscanVersion)
                            generator.writeStringField("interpro-version", dbReleases?.interpro?.version)
                            generator.writeFieldName("results")
                            generator.writeStartArray()
                            writeProtein(proteinMd5, proteinMatches, generator, seqData)
                            generator.writeEndArray()
                            generator.writeEndObject()
                            generator.writeRaw('\n')
                        }
                    }
                }
            } else {
                DefaultPrettyPrinter pp = new DefaultPrettyPrinter()
                pp.indentArraysWith(DefaultIndenter.SYSTEM_LINEFEED_INSTANCE)
                generator.setPrettyPrinter(pp)
                generator.writeStartObject()

                generator.writeStringField("interproscan-version", iprscanVersion)
                generator.writeStringField("interpro-version", dbReleases?.interpro?.version)
                generator.writeFieldName("results")
                generator.writeStartArray()  // start of results [...
                Set<String> seenNucleicMd5s = new HashSet<>()
                inputPaths.each { inputPath ->
                    Map proteins = mapper.readValue(new File(inputPath), Map)
                    if (isNucleic) {  // input was nucleic acid sequence
                        def (nucleicToProteinMd5, ntSeqDataMap, orfDataMap) = 
                            db.retrieveAllNucleicSequenceData(proteins.keySet() as List)
                        nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                            if (!seenNucleicMd5s.contains(nucleicMd5)) {
                                def seqData = ntSeqDataMap[nucleicMd5]
                                def protMd5ToOrfs = orfDataMap[nucleicMd5]
                                writeNucleic(nucleicMd5, proteinMd5s, proteins, generator, seqData, protMd5ToOrfs)
                                seenNucleicMd5s.add(nucleicMd5)
                            }
                        }
                    } else {  // input was protein sequences
                        def allMd5s = proteins.keySet() as List
                        def md5ToSeqData = [:]
                        allMd5s.collate(1000).each { batch ->
                            def result = db.proteinMd5ToProteinSeqs(batch)
                            md5ToSeqData.putAll(result)
                        }

                        proteins.each { String proteinMd5, Map proteinMatches ->
                            def seqData = md5ToSeqData[proteinMd5]
                            writeProtein(proteinMd5, proteinMatches, generator, seqData)
                        }
                    }
                }
                generator.writeEndArray() // end of "results" ...]
                generator.writeEndObject()
            }
        }

    }

    static def streamJson(String filePath, Closure closure) {
        try (FileWriter fileWriter = new FileWriter(new File(filePath))) {
            ObjectMapper mapper = new ObjectMapper()
            try (JsonGenerator generator = mapper.getFactory().createGenerator(fileWriter)) {
                closure.call(generator)
            }    
        }
    }

    static void writeNucleic(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, List seqData, Map protMd5ToOrfs) {
        /* Write data for an input nucleic acid sequence, and then the matches for its associated ORFs
        {"sequence: nt seq, "md5": nt md5,
        "crossReferences": [{ntSeqData}, {ntSeqData}],
        "openReadingFrames": [{protein}, {protein}, {protein}]}
        There may be multiple nt seq Ids associated with the same nt seq, use the first entry to get the seq. */
        jsonWriter.writeStartObject()

        // 1. {"sequence": seq, "md5": ntMd5}
        String sequence = seqData[0].sequence
        jsonWriter.writeStringField("sequence", sequence)
        jsonWriter.writeStringField("md5", nucleicMd5)

        // 2. {..., "crossReferences": [{ntSeqXref}, {ntSeqXref}]}
        jsonWriter.writeFieldName("crossReferences")
        writeXref(seqData, jsonWriter)

        // 3. {..., "openReadingFrames": [{protein}, {protein}]}
        jsonWriter.writeFieldName("openReadingFrames")
        writeOpenReadingFrames(nucleicMd5, proteinMd5s, proteinMatches, jsonWriter, protMd5ToOrfs)

        jsonWriter.writeEndObject()
    }

    static void writeOpenReadingFrames(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, Map protMd5ToOrfs){
        def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

        jsonWriter.writeStartArray()
        proteinMd5s.each { String proteinMd5 ->
            // a proteinSeq/Md5 may be associated with multiple nt md5s/seq, only pull the data where the nt md5/seq is relevant
            def proteinSeqData = protMd5ToOrfs[proteinMd5]
            proteinSeqData.each { row ->
                def proteinSource = SOURCE_NT_PATTERN.matcher(row.description)
                assert proteinSource.matches()
                jsonWriter.writeStartObject()
                jsonWriter.writeNumberField("start", proteinSource.group(1) as int)
                jsonWriter.writeNumberField("end", proteinSource.group(2) as int)
                jsonWriter.writeStringField("strand", (proteinSource.group(3) as int) < 4 ? "SENSE" : "ANTISENSE")
                jsonWriter.writeFieldName("protein")
                writeProtein(proteinMd5, proteinMatches[proteinMd5], jsonWriter, proteinSeqData)
                jsonWriter.writeEndObject()
            }
        }
        jsonWriter.writeEndArray()
    }

    static void writeProtein(String proteinMd5, Map proteinMatches, JsonGenerator jsonWriter, List proteinSeqData) {
        /* Write data for a query protein sequence and its matches:
        { "sequence": sequence, "md5": proteinMd5, "matches": [], "xrefs": []}
        There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */
        jsonWriter.writeStartObject()

        // 1. {"sequence": seq, "md5": proteinMd5}
        String sequence = proteinSeqData[0].sequence
        jsonWriter.writeStringField("sequence", sequence)
        jsonWriter.writeStringField("md5", proteinMd5)

        // 2. {..., "matches": [{match}, {match}, {match}]}
        jsonWriter.writeFieldName("matches")
        jsonWriter.writeStartArray()
        proteinMatches.each { String modelAcc, Map match->
            writeMatch(proteinMd5, match, jsonWriter)
        }
        jsonWriter.writeEndArray()

        // 3. {..., "xref": [{xref}, {xref}, {xref}]}
        jsonWriter.writeFieldName("xref")
        writeXref(proteinSeqData, jsonWriter)
        jsonWriter.writeEndObject()
    }

    static void writeMatch(String proteinMd5, Map match, JsonGenerator jsonWriter) {
        // Write out an individual match to an array of matches. The structure is dependent on the application.
        String appl = (match.source == "InterPro-N" ? "InterPro-N" 
                                                    : match.signature.signatureLibraryRelease.library
                                                    ).toLowerCase()
        switch (appl) {
            case "antifam":
                writeDefault(match, jsonWriter)
                break
            case "cath-gene3d":
                writeDefault(match, jsonWriter)
                break
            case "cath-funfam":
                writeDefault(match, jsonWriter)
                break
            case "cdd":
                writeCDD(match, jsonWriter)
                break
            case "coils":
                writeMinimalist(match, jsonWriter)
                break
            case "hamap":
                writeHAMAP(match, jsonWriter)
                break
            case "interpro-n":
                writeInterProN(match, jsonWriter)
                break
            case "mobidb-lite":
                writeMobiDBlite(match, jsonWriter)
                break
            case "ncbifam":
                writeDefault(match, jsonWriter)
                break
            case "panther":
                writePANTHER(match, jsonWriter)
                break
            case "pfam":
                writeDefault(match, jsonWriter)
                break
            case "phobius":
                writeMinimalist(match, jsonWriter)
                break
            case "pirsf":
                writeDefault(match, jsonWriter)
                break
            case "pirsr":
                writePIRSR(match, jsonWriter)
                break
            case "prints":
                writePRINTS(match, jsonWriter)
                break
            case "prosite patterns":
                writePROSITEpatterns(match, jsonWriter)
                break
            case "prosite profiles":
                writePROSITEprofiles(match, jsonWriter)
                break
            case "sfld":
                writeSFLD(match, jsonWriter)
                break
            case "signalp":
                writeSignalp(match, jsonWriter)
                break
            case "smart":
                writeSMART(match, jsonWriter)
                break
            case "superfamily":
                writeSUPERFAMILY(match, jsonWriter)
                break
            case "deeptmhmm":
                writeMinimalist(match, jsonWriter)
                break
            case "tmbed":
                writeMinimalist(match, jsonWriter)
                break
            default:
                throw new UnsupportedOperationException("Unknown application '${appl}' for query protein with MD5 ${proteinMd5}")
        }
    }

    static void writeDefault(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "score"    : match.score,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "hmmBounds"         : Location.getHmmBounds(loc.hmmBounds),
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "envelopeStart"     : loc.envelopeStart,
                    "envelopeEnd"       : loc.envelopeEnd,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeDefaultNoHmmBounds(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "score"    : match.score,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "envelopeStart"     : loc.envelopeStart,
                    "envelopeEnd"       : loc.envelopeEnd,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeCDD(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "location-fragments": formatFragments(loc.fragments),
                    "sites"             : loc.sites.collect { site ->
                        [
                            "description"  : site.description,
                            "numLocations" : site.numLocations,
                            "siteLocations": site.siteLocations
                        ]
                    }
                ]
            }
        ])
    }

    static void writeMinimalist(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeHAMAP(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "score"             : loc.score,
                    "cigarAlignment"    : loc.cigarAlignment,
                    "alignment"         : loc.targetAlignment,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeMobiDBlite(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "location-fragments": formatFragments(loc.fragments),
                    "sequence-feature"  : loc.sequenceFeature
                ]
            }
        ])
    }

    static void writePANTHER(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature"      : match.signature,
            "model-ac"       : match.treegrafter.subfamilyAccession ?: match.modelAccession,
            "name"           : match.treegrafter.subfamilyDescription,
            "evalue"         : match.evalue,
            "score"          : match.score,
            "proteinClass"   : match.treegrafter.proteinClass,
            "graftPoint"     : match.treegrafter.graftPoint,
            "ancestralNode"  : match.treegrafter.ancestralNodeID,
            "goXRefs"        : match.treegrafter.goXRefs,
            "source"         : match.source,
            "locations"      : match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "hmmBounds"         : Location.getHmmBounds(loc.hmmBounds),
                    "envelopeStart"     : loc.envelopeStart,
                    "envelopeEnd"       : loc.envelopeEnd,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writePIRSR(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "score"    : match.score,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "score"             : loc.score,
                    "envelopeStart"     : loc.envelopeStart,
                    "envelopeEnd"       : loc.envelopeEnd,
                    "location-fragments": formatFragments(loc.fragments),
                    "sites"             : loc.sites.collect { site ->
                        [
                            "description"  : site.description,
                            "numLocations" : site.numLocations,
                            "siteLocations": site.siteLocations.collect { siteLoc ->
                                [
                                    "start"  : siteLoc.start,
                                    "end"    : siteLoc.end,
                                    "residue": siteLoc.residue
                                ]
                            },
                        ]
                    } // end of "sites"
                ]
            } // end of "locations"
        ])
    }

    static void writePRINTS(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "graphscan": match.graphscan,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "pvalue"            : loc.pvalue,
                    "score"             : loc.score,
                    "motifNumber"       : loc.motifNumber,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writePROSITEpatterns(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "level"             : loc.level,
                    "cigarAlignment"    : loc.cigarAlignment,
                    "alignment"         : loc.targetAlignment,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writePROSITEprofiles(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "score"             : loc.score,
                    "cigarAlignment"    : loc.cigarAlignment,
                    "alignment"         : loc.targetAlignment,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeSignalp(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "score"             : loc.score,
                    "location-fragments": formatFragments(loc.fragments),
                ]
            }
        ])
    }

    static void writeSFLD(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "score"    : match.score,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "envelopeStart"     : loc.envelopeStart,
                    "envelopeEnd"       : loc.envelopeEnd,
                    "location-fragments": formatFragments(loc.fragments),
                    "sites"             : loc.sites.collect { site ->
                        [
                            "description": site.description,
                            "numLocations": site.numLocations,
                            "siteLocations": site.siteLocations.collect { siteLoc ->
                                [
                                    "start"  : siteLoc.start,
                                    "end"    : siteLoc.end,
                                    "residue": siteLoc.residue
                                ]
                            },
                        ]
                    } // end of "sites"
                ]
            } // end of "locations"
        ])
    }

    static void writeSMART(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "evalue"   : match.evalue,
            "score"    : match.score,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmStart"          : loc.hmmStart,
                    "hmmEnd"            : loc.hmmEnd,
                    "hmmLength"         : loc.hmmLength,
                    "hmmBounds"         : Location.getHmmBounds(loc.hmmBounds),
                    "evalue"            : loc.evalue,
                    "score"             : loc.score,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeSUPERFAMILY(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "representative"    : loc.representative,
                    "hmmLength"         : loc.hmmLength,
                    "evalue"            : loc.evalue,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static void writeInterProN(Map match, JsonGenerator jsonWriter) {
        jsonWriter.writeObject([
            "signature": match.signature,
            "model-ac" : match.modelAccession,
            "source"   : match.source,
            "locations": match.locations.collect { loc ->
                [
                    "start"             : loc.start,
                    "end"               : loc.end,
                    "score"             : loc.score,
                    "representative"    : loc.representative,
                    "location-fragments": formatFragments(loc.fragments)
                ]
            }
        ])
    }

    static List<Map> formatFragments(fragments) {
        return fragments.collect { frag ->
            [
                "start"    : frag.start,
                "end"      : frag.end,
                "dc-status": frag.dcStatus
            ]
        }
    }

    static void writeXref(List seqData, JsonGenerator jsonWriter) {
        /* "xref"/"crossReferences" : [ {
            "name" : "tr|A0A011PH51|A0A011PH51_9PROT OX=1454000",
            "id" : "tr|A0A011PH51|A0A011PH51_9PROT"
        } ] */
        jsonWriter.writeStartArray()
        seqData.each { row ->
            // jsonWrite.writeObject([name: "$seqId $seqDesc"]) does not correctly handle the formatted str
            jsonWriter.writeStartObject()
            jsonWriter.writeStringField("name", "${row.id} ${row.description}".trim())
            jsonWriter.writeStringField("id", row.id)
            jsonWriter.writeEndObject()
        }
        jsonWriter.writeEndArray()
    }
}
