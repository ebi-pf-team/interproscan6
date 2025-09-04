import com.fasterxml.jackson.core.JsonFactory
import com.fasterxml.jackson.core.JsonGenerator
import com.fasterxml.jackson.core.io.SerializedString
import com.fasterxml.jackson.core.util.DefaultPrettyPrinter
import com.fasterxml.jackson.core.util.DefaultIndenter
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.SerializationFeature
import groovy.json.JsonException
import java.util.regex.Pattern

import Match

process WRITE_JSON {
    label    'tiny'
    executor 'local'

    input:
    val matches_files  // {query prot seq md5: {model acc: match}}
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val db_releases
    val jsonlines

    exec:
    ObjectMapper jacksonMapper = new ObjectMapper()
    SeqDB db = new SeqDB(seq_db_file.toString())

    streamJson(output_file.toString(), jacksonMapper, jsonlines) { JsonGenerator generator ->
        if (jsonlines) {
            generator.setRootValueSeparator(new SerializedString(''))
            Set<String> seenNucleicMd5s = new HashSet<>()
            matches_files.each { matchFile ->
                Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)

                if (nucleic) {
                    nucleicToProteinMd5 = db.groupProteins(proteins)
                    nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                        generator.writeStartObject()
                        generator.writeStringField("interproscan-version", interproscan_version)
                        generator.writeStringField("interpro-version", db_releases?.interpro?.version)
                        generator.writeFieldName("results")
                        generator.writeStartArray()
                        if (!seenNucleicMd5s.contains(nucleicMd5)) {
                            writeNucleic(nucleicMd5, proteinMd5s, proteins, generator, db)
                            seenNucleicMd5s.add(nucleicMd5)
                        }
                        generator.writeEndArray()
                        generator.writeEndObject()
                        generator.writeRaw('\n')
                    }
                } else {
                    proteins.each { String proteinMd5, Map proteinMatches ->
                        generator.writeStartObject()
                        generator.writeStringField("interproscan-version", interproscan_version)
                        generator.writeStringField("interpro-version", db_releases?.interpro?.version)
                        generator.writeFieldName("results")
                        generator.writeStartArray()
                        writeProtein(proteinMd5, proteinMatches, generator, db)
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

            generator.writeStringField("interproscan-version", interproscan_version)
            generator.writeStringField("interpro-version", db_releases?.interpro?.version)
            generator.writeFieldName("results")
            generator.writeStartArray()  // start of results [...
            Set<String> seenNucleicMd5s = new HashSet<>()
            matches_files.each { matchFile ->
                Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
                if (nucleic) {  // input was nucleic acid sequence
                    nucleicToProteinMd5 = db.groupProteins(proteins)
                    nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                        if (!seenNucleicMd5s.contains(nucleicMd5)) {
                            writeNucleic(nucleicMd5, proteinMd5s, proteins, generator, db)
                            seenNucleicMd5s.add(nucleicMd5)
                        }
                    }
                } else {  // input was protein sequences
                    proteins.each { String proteinMd5, Map proteinMatches ->
                        writeProtein(proteinMd5, proteinMatches, generator, db)
                    }
                }
            }
            generator.writeEndArray() // end of "results" ...]
            generator.writeEndObject()
        }
    }
}

def writeNucleic(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, SeqDB db) {
    /* Write data for an input nucleic acid sequence, and then the matches for its associated ORFs
    {"sequence: nt seq, "md5": nt md5,
    "crossReferences": [{ntSeqData}, {ntSeqData}],
    "openReadingFrames": [{protein}, {protein}, {protein}]}
    There may be multiple nt seq Ids associated with the same nt seq, use the first entry to get the seq. */
    jsonWriter.writeStartObject()

    // 1. {"sequence": seq, "md5": ntMd5}
    ntSeqData = db.nucleicMd5ToNucleicSeq(nucleicMd5)
    String sequence = ntSeqData[0].sequence
    jsonWriter.writeStringField("sequence", sequence)
    jsonWriter.writeStringField("md5", nucleicMd5)

    // 2. {..., "crossReferences": [{ntSeqXref}, {ntSeqXref}]}
    jsonWriter.writeFieldName("crossReferences")
    writeXref(ntSeqData, jsonWriter)

    // 3. {..., "openReadingFrames": [{protein}, {protein}]}
    jsonWriter.writeFieldName("openReadingFrames")
    writeOpenReadingFrames(nucleicMd5, proteinMd5s, proteinMatches, jsonWriter, db)

    jsonWriter.writeEndObject()
}

def writeOpenReadingFrames(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, JsonGenerator jsonWriter, SeqDB db){
    def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

    jsonWriter.writeStartArray()
    proteinMd5s.each { String proteinMd5 ->
        // a proteinSeq/Md5 may be associated with multiple nt md5s/seq, only pull the data where the nt md5/seq is relevant
        proteinSeqData = db.getOrfSeq(proteinMd5, nucleicMd5)
        proteinSeqData.each { row ->
            def proteinSource = SOURCE_NT_PATTERN.matcher(row.description)
            assert proteinSource.matches()
            jsonWriter.writeStartObject()
            jsonWriter.writeNumberField("start", proteinSource.group(1) as int)
            jsonWriter.writeNumberField("end", proteinSource.group(2) as int)
            jsonWriter.writeStringField("strand", (proteinSource.group(3) as int) < 4 ? "SENSE" : "ANTISENSE")
            jsonWriter.writeFieldName("protein")
            writeProtein(proteinMd5, proteinMatches[proteinMd5], jsonWriter, db)
            jsonWriter.writeEndObject()
        }
    }
    jsonWriter.writeEndArray()
}

def writeProtein(String proteinMd5, Map proteinMatches, JsonGenerator jsonWriter, SeqDB db) {
    /* Write data for a query protein sequence and its matches:
    { "sequence": sequence, "md5": proteinMd5, "matches": [], "xrefs": []}
    There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */
    jsonWriter.writeStartObject()

    // 1. {"sequence": seq, "md5": proteinMd5}
    proteinSeqData = db.proteinMd5ToProteinSeq(proteinMd5)
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

def writeMatch(String proteinMd5, Map match, JsonGenerator jsonWriter) {
    // Write out an individual match to an array of matches. The structure is dependent on the memberDB.
    String memberDB = match.signature.signatureLibraryRelease.library.toLowerCase() ?: ""
    switch (memberDB) {
        case "antifam":
            writeDefault(match, jsonWriter)
            break
        case "cath-gene3d":
            writeDefault(match, jsonWriter)
            break
        case "cath-funfam":
        case "funfam":  // use groovy case fall to allow multiple options
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
        case "mobidb lite":
        case "mobidb-lite":
        case "mobidb_lite":  // use groovy case fall to allow multiple options
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
            writePirsr(match, jsonWriter)
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
        case "tmhmm":
        case "deeptmhmm":
            writeMinimalist(match, jsonWriter)
            break
        case "tmbed":
            writeMinimalist(match, jsonWriter)
            break
        default:
            throw new UnsupportedOperationException("Unknown database '${memberDB}' for query protein with MD5 ${proteinMd5}")
    }
}

def writeDefault(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"    : match.score,
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

def writeDefaultNoHmmBounds(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"    : match.score,
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

def writeCDD(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"    : match.score,
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

def writeMinimalist(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
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

def writeHAMAP(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
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

def writeMobiDBlite(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
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

def writePANTHER(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature"      : match.signature,
        "model-ac"       : match.treegrafter.subfamilyAccession ?: match.modelAccession,
        "name"           : match.treegrafter.subfamilyDescription,
        "evalue"         : match.evalue,
        "score"          : match.score,
        "proteinClass"   : match.treegrafter.proteinClass,
        "graftPoint"     : match.treegrafter.graftPoint,
        "ancestralNode": match.treegrafter.ancestralNodeID,
        "goXRefs"        : match.treegrafter.goXRefs,
        "locations"      : match.locations.collect { loc ->
            [
                "start"             : loc.start,
                "end"               : loc.end,
                "representative"    : loc.representative,
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

def writePirsr(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"    : match.score,
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

def writePRINTS(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "graphscan": match.graphscan,
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

def writePROSITEpatterns(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
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

def writePROSITEprofiles(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
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

def writeSignalp(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "locations": match.locations.collect { loc ->
            [
                "start"             : loc.start,
                "end"               : loc.end,
                "representative"    : loc.representative,
                "pvalue"            : loc.score,
                "location-fragments": formatFragments(loc.fragments),
            ]
        }
    ])
}

def writeSFLD(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"    : match.score,
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

def writeSMART(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "score"     : match.score,
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

def writeSUPERFAMILY(Map match, JsonGenerator jsonWriter) {
    jsonWriter.writeObject([
        "signature": match.signature,
        "model-ac" : match.modelAccession,
        "evalue"   : match.evalue,
        "locations": match.locations.collect { loc ->
            [
                "start"             : loc.start,
                "end"               : loc.end,
                "representative"    : loc.representative,
                "hmmLength"         : loc.hmmLength,
                "location-fragments": formatFragments(loc.fragments)
            ]
        }
    ])
}

def formatFragments(fragments) {
    return fragments.collect { frag ->
        [
            "start"    : frag.start,
            "end"      : frag.end,
            "dc-status": frag.dcStatus
        ]
    }
}

def writeXref(seqData, JsonGenerator jsonWriter) {
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

def streamJson(String filePath, ObjectMapper mapper, boolean jsonlines, Closure closure) {
    FileWriter fileWriter = null
    JsonGenerator generator = null
    try {
        JsonFactory factory = mapper.getFactory()
        fileWriter = new FileWriter(new File(filePath))
        generator = factory.createGenerator(fileWriter)
        closure.call(generator)  // Call the closure to write key-value pairs
    } catch (IOException e) {
        throw new JsonException("IO error writing file: $filePath\nException: $e\nCause: ${e.getCause()}", e)
    } catch (Exception e) {
        throw new Exception("Error occurred when writing Json file $filePath\nException: $e\nCause: ${e.getCause()}", e)
    } finally {
        if (generator != null) {
            generator.close()
        }
        if (fileWriter != null) {
            fileWriter.close()
        }
    }
}