import groovy.xml.StreamingMarkupBuilder
import java.io.FileWriter
import java.util.regex.Pattern
import javax.xml.stream.XMLOutputFactory
import javax.xml.stream.XMLStreamWriter
import com.fasterxml.jackson.databind.ObjectMapper

import java.time.format.DateTimeFormatter
import java.time.LocalDate

import Match

process WRITE_XML {
    label    'tiny'
    executor 'local'

    input:
    val matches_files  // {query prot seq md5: {model acc: match}}
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val db_releases

    exec:
    def db = new SeqDBQuery(seq_db_file.toString())

    // Use BufferedWriter with custom buffer size
    def bufferedWriter = new BufferedWriter(new FileWriter(output_file), 1024 * 1024) // 1MB buffer
    
    def timingFile = new File("${output_file}_timing.log")
    timingFile.text = "timestamp,operation,duration_ms,details\n"

    try {
        def factory = XMLOutputFactory.newInstance()
        def writer = factory.createXMLStreamWriter(bufferedWriter)
        
        writer.writeStartDocument("UTF-8", "1.0")
        writer.writeStartElement("results")
        writer.writeAttribute("interproscan-version", interproscan_version)
        writer.writeAttribute("interpro-version", db_releases?.interpro?.version ?: "")

        Set<String> seenNucleicMd5s = new HashSet<>()
        def proteinCount = 0
        def BUFFER_SIZE = 5000  // Flush every 5000 proteins

        matches_files.each { matchFile ->
            def jsonLoadStartTime = System.currentTimeMillis()
            Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
            def jsonLoadEndTime = System.currentTimeMillis()
            def jsonLoadDuration = jsonLoadEndTime - jsonLoadStartTime

            def timestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
            def fileSize = new File(matchFile.toString()).length()
            timingFile.append("${timestamp},json_load,${jsonLoadDuration},\"${new File(matchFile.toString()).name} (${fileSize} bytes, ${proteins.size()} proteins)\"\n")
            
            if (nucleic) {
                proteinCount = processNucleotidesBulkBuffered(db, proteins, seenNucleicMd5s, writer, bufferedWriter, BUFFER_SIZE, proteinCount, timingFile)
            } else {
                def proteinQueryStartTime = System.currentTimeMillis()
                def proteinMd5List = proteins.keySet().toList()
                Map<String, List> seqData = db.proteinMd5sToProteinSeqs(proteinMd5List)
                def proteinQueryEndTime = System.currentTimeMillis()
                def proteinQueryDuration = proteinQueryEndTime - proteinQueryStartTime
                
                def proteinQueryTimestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
                timingFile.append("${proteinQueryTimestamp},db_query,${proteinQueryDuration},\"${proteinMd5List.size()} proteins from ${new File(matchFile.toString()).name}\"\n")
                
                def proteinProcessingStartTime = System.currentTimeMillis()
                
                for (entry in proteins) {
                    String proteinMd5 = entry.key
                    Map proteinMatches = entry.value
                    List proteinSeqData = seqData[proteinMd5]
                    
                    def individualProteinStartTime = System.currentTimeMillis()
                    addProteinNodesDirect(proteinMd5, proteinMatches, proteinSeqData, writer)
                    def individualProteinEndTime = System.currentTimeMillis()
                    def individualProteinDuration = individualProteinEndTime - individualProteinStartTime
                    
                    proteinCount++
                    
                    if (proteinCount % BUFFER_SIZE == 0) {
                        def flushStartTime = System.currentTimeMillis()
                        writer.flush()
                        bufferedWriter.flush()
                        def flushEndTime = System.currentTimeMillis()
                        def flushDuration = flushEndTime - flushStartTime
                        
                        def flushTimestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
                        timingFile.append("${flushTimestamp},buffer_flush,${flushDuration},\"Flushed after ${proteinCount} proteins\"\n")
                        println "Flushed buffer after processing ${proteinCount} proteins..."
                    }
                    
                    def proteinTimestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
                    def matchCount = proteinMatches.size()
                    timingFile.append("${proteinTimestamp},protein_processing,${individualProteinDuration},\"${proteinMd5} with ${matchCount} matches\"\n")
                }
                
                def proteinLoopEndTime = System.currentTimeMillis()
                def proteinLoopDuration = proteinLoopEndTime - proteinProcessingStartTime
                
                def proteinLoopTimestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
                timingFile.append("${proteinLoopTimestamp},proteins_loop_total,${proteinLoopDuration},\"${proteins.size()} proteins from ${new File(matchFile.toString()).name}\"\n")
            }
        }

        writer.writeEndElement() // results
        writer.writeEndDocument()
        
        // Final flush
        writer.flush()
        bufferedWriter.flush()
        writer.close()

    } finally {
        bufferedWriter.close()
        db.close()
    }
}

def processNucleotidesBulkBuffered(SeqDBQuery db, Map proteins, Set seenNucleicMd5s, XMLStreamWriter writer, 
                                   BufferedWriter bufferedWriter, int bufferSize, int currentCount, File timingFile) {
    Set<String> allProteinMd5s = proteins.keySet().toSet()
    Map<String, Set<String>> nucleicToProteinMd5 = db.groupProteinsBulk(allProteinMd5s)
    
    List<String> newNucleicMd5s = nucleicToProteinMd5.keySet().findAll { !seenNucleicMd5s.contains(it) }
    seenNucleicMd5s.addAll(newNucleicMd5s)
    
    if (newNucleicMd5s.isEmpty()) {
        return currentCount
    }

    Map<String, List> nucleotideSeqData = db.nucleicMd5sToNucleicSeqs(newNucleicMd5s)

    List<String> relevantProteinMd5s = []
    newNucleicMd5s.each { nucleicMd5 ->
        relevantProteinMd5s.addAll(nucleicToProteinMd5[nucleicMd5])
    }
    Map<String, Map<String, List>> orfSeqData = db.getOrfSeqsBulk(relevantProteinMd5s, newNucleicMd5s)

    def processedCount = currentCount
    newNucleicMd5s.each { String nucleicMd5 ->
        def ntSeqData = nucleotideSeqData[nucleicMd5]
        if (!ntSeqData) return

        Set<String> proteinMd5sForNucleic = nucleicToProteinMd5[nucleicMd5]

        addNucleotideNodesDirect(nucleicMd5, proteinMd5sForNucleic, proteins, ntSeqData, orfSeqData, writer)
        processedCount++
        
        if (processedCount % bufferSize == 0) {
            def flushStartTime = System.currentTimeMillis()
            writer.flush()
            bufferedWriter.flush()
            def flushEndTime = System.currentTimeMillis()
            def flushDuration = flushEndTime - flushStartTime
            
            def flushTimestamp = new Date().format("yyyy-MM-dd HH:mm:ss.SSS")
            timingFile.append("${flushTimestamp},buffer_flush,${flushDuration},\"Flushed after ${processedCount} nucleotides\"\n")
            println "Flushed buffer after processing ${processedCount} nucleotides..."
        }
    }
    
    return processedCount
}

def addNucleotideNodesDirect(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, List ntSeqData, 
                           Map orfSeqData, XMLStreamWriter writer) {
    def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)
    
    String sequence = ntSeqData[0].sequence

    writer.writeStartElement("nucleotide")
    
    // <sequence md5="" sequence </sequence>
    writer.writeStartElement("sequence")
    writer.writeAttribute("md5", nucleicMd5)
    writer.writeCharacters(sequence)
    writer.writeEndElement()

    // <xref id="id" name="id desc"/>
    writeXrefDirect(ntSeqData, writer)

    // <orf end="", start="", strand="">
    proteinMd5s.each { proteinMd5 ->
        def proteinSeqData = orfSeqData[proteinMd5]?.get(nucleicMd5)
        if (proteinSeqData) {
            proteinSeqData.each { row ->
                def proteinSource = SOURCE_NT_PATTERN.matcher(row.description)
                if (proteinSource.matches()) {
                    writer.writeStartElement("orf")
                    writer.writeAttribute("start", proteinSource.group(1))
                    writer.writeAttribute("end", proteinSource.group(2))
                    writer.writeAttribute("strand", proteinSource.group(3) as int < 4 ? "SENSE" : "ANTISENSE")
                    
                    addProteinNodesDirect(proteinMd5, proteinMatches[proteinMd5], [row], writer)
                    
                    writer.writeEndElement() // orf
                }
            }
        }
    }
    
    writer.writeEndElement() // nucleotide-sequence
}

def addProteinNodesDirect(String proteinMd5, Map proteinMatches, List proteinSeqData, XMLStreamWriter writer) {
    if (!proteinSeqData || proteinSeqData.size() == 0) return
    
    String sequence = proteinSeqData[0].sequence
    
    writer.writeStartElement("protein")
    
    // <sequence md5="" sequence </sequence>
    writer.writeStartElement("sequence")
    writer.writeAttribute("md5", proteinMd5)
    writer.writeCharacters(sequence)
    writer.writeEndElement()

    // <xref id="id", name="id desc"/>
    writeXrefDirect(proteinSeqData, writer)

    // <matches>
    writer.writeStartElement("matches")
    for (entry in proteinMatches) {
        addMatchNodeDirect(proteinMd5, entry.value, writer)
    }
    writer.writeEndElement() // matches
    
    writer.writeEndElement() // protein
}

def addMatchNodeDirect(String proteinMd5, Map match, XMLStreamWriter writer) {
    String memberDB = match.signature.signatureLibraryRelease.library.toLowerCase() ?: ""

    def matchNodeAttributes = getMatchNodeAttributes(memberDB, match)

    writer.writeStartElement("match")
    
    // Write match attributes
    if (matchNodeAttributes) {
        matchNodeAttributes.each { key, value ->
            if (value != null) {
                writer.writeAttribute(key.toString(), value.toString())
            }
        }
    }
    
    // <signature>
    def signatureNodeAttributes = fmtSignatureNode(match)
    writer.writeStartElement("signature")
    signatureNodeAttributes.each { key, value ->
        if (value != null) {
            writer.writeAttribute(key.toString(), value.toString())
        }
    }
    
    // <entry> if exists
    if (match.signature.entry) {
        addEntryNodeDirect(match.signature.entry, writer)
    }
    
    // <signature-library-release>
    writer.writeStartElement("signature-library-release")
    writer.writeAttribute("library", match.signature.signatureLibraryRelease.library)
    writer.writeAttribute("version", match.signature.signatureLibraryRelease.version)
    writer.writeEndElement()
    
    writer.writeEndElement() // signature
    
    // <model-ac>
    writer.writeStartElement("model-ac")
    if (memberDB == "panther") {
        writer.writeCharacters(match.treegrafter?.subfamilyAccession ?: match.modelAccession)
        writer.writeEndElement() // model-ac
        
        // PANTHER go-xrefs
        match.treegrafter?.goXRefs?.each { goXref ->
            writer.writeStartElement("go-xref")
            writer.writeAttribute("category", goXref.category)
            writer.writeAttribute("db", goXref.databaseName)
            writer.writeAttribute("id", goXref.id)
            writer.writeAttribute("name", goXref.name)
            writer.writeEndElement()
        }
    } else {
        writer.writeCharacters(match.modelAccession)
        writer.writeEndElement() // model-ac
    }
    
    // <locations>
    addLocationNodesDirect(memberDB, proteinMd5, match, writer)
    
    writer.writeEndElement() // match
}

def getMatchNodeAttributes(String memberDB, Map match) {
    switch (memberDB) {
        case "antifam":
        case "cath-gene3d":
        case "cath-funfam":
        case "funfam":
        case "cdd":
        case "hamap":
        case "ncbifam":
        case "pfam":
        case "phobius":
        case "pirsf":
        case "pirsr":
        case "prosite patterns":
        case "prosite profiles":
        case "sfld":
        case "smart":
            return fmtDefaultMatchNode(match)
        case "panther":
            return fmtPantherMatchNode(match)
        case "prints":
            return fmtPrintsMatchNode(match)
        case "superfamily":
            return fmtSuperfamilyMatchNode(match)
        case "coils":
        case "mobidb lite":
        case "mobidb-lite":
        case "mobidb_lite":
        case "signalp":
        case "tmhmm":
        case "deeptmhmm":
        case "tmbed":
            return null
        default:
            throw new UnsupportedOperationException("Unknown database '${memberDB}' for query protein with MD5 ${proteinMd5}")
    }
}

def addEntryNodeDirect(Map entry, XMLStreamWriter writer) {
    writer.writeStartElement("entry")
    writer.writeAttribute("ac", entry.accession)
    writer.writeAttribute("desc", entry.description)
    writer.writeAttribute("name", entry.name)
    writer.writeAttribute("type", entry.type)
    
    // GO cross-references
    entry.goXRefs?.each { goXref ->
        writer.writeStartElement("go-xref")
        writer.writeAttribute("category", goXref.category)
        writer.writeAttribute("db", goXref.databaseName)
        writer.writeAttribute("id", goXref.id)
        writer.writeAttribute("name", goXref.name)
        writer.writeEndElement()
    }
    
    // Pathway cross-references
    entry.pathwayXRefs?.each { pathwayXref ->
        writer.writeStartElement("pathway-xref")
        writer.writeAttribute("db", pathwayXref.databaseName)
        writer.writeAttribute("id", pathwayXref.id)
        writer.writeAttribute("name", pathwayXref.name)
        writer.writeEndElement()
    }
    
    writer.writeEndElement() // entry
}

def addLocationNodesDirect(String memberDB, String proteinMd5, Map match, XMLStreamWriter writer) {
    writer.writeStartElement("locations")
    
    match.locations.each { loc ->
        def locationAttributes = getLocationAttributes(memberDB, match, loc)
        
        writer.writeStartElement("location")
        locationAttributes.each { key, value ->
            if (value != null) {
                writer.writeAttribute(key.toString(), value.toString())
            }
        }
        
        // fragments if they exist
        if (loc.containsKey("fragments") && loc.fragments?.size() > 0) {
            writer.writeStartElement("location-fragments")
            loc.fragments.each { frag ->
                writer.writeStartElement("fragment")
                writer.writeAttribute("start", frag.start.toString())
                writer.writeAttribute("end", frag.end.toString())
                writer.writeAttribute("dc-status", frag.dcStatus)
                writer.writeEndElement()
            }
            writer.writeEndElement() // location-fragments
        }
        
        // alignment info for specific databases
        if (memberDB in ["hamap", "prosite patterns", "prosite profiles"]) {
            writer.writeStartElement("alignment")
            writer.writeCharacters(loc.targetAlignment ?: "")
            writer.writeEndElement()
            
            writer.writeStartElement("cigar-alignment")
            writer.writeCharacters(loc.cigarAlignment ?: "")
            writer.writeEndElement()
        }
        
        // sites if they exist
        if (loc.containsKey("sites") && loc.sites?.size() > 0) {
            addSiteNodesDirect(loc.sites, writer)
        }
        
        writer.writeEndElement() // location
    }
    
    writer.writeEndElement() // locations
}

def getLocationAttributes(String memberDB, Map match, Map loc) {
    switch (memberDB) {
        case "antifam":
        case "cath-funfam":
        case "funfam":
        case "cath-gene3d":
        case "gene3d":
        case "ncbifam":
        case "pfam":
        case "pirsf":
            return fmtDefaultLocationNode(loc)
        case "cdd":
            return fmtCddLocationNode(match, loc)
        case "coils":
        case "tmhmm":
        case "deeptmhmm":
        case "tmbed":
            return fmMinimalistLoctationNode(loc)
        case "hamap":
        case "phobius":
        case "prosite profiles":
            return fmtMinimalistLocationNode(loc)
        case "mobidb lite":
        case "mobidb-lite":
        case "mobidb_lite":
            return fmtMobidbLiteLocationNode(loc)
        case "panther":
            return fmtPantherLocationNode(loc)
        case "prints":
            return fmtPrintsLocationNode(loc)
        case "prosite patterns":
            return fmtPrositePatternsLocationNode(loc)
        case "pirsr":
        case "sfld":
            return fmtDefaultNoHbLocationNode(loc)
        case "signalp":
            return fmtSignalpLocationNode(loc)
        case "smart":
            return fmtSmartLocationNode(loc)
        case "superfamily":
            return fmtSuperfamilyLocationNode(loc)
        default:
            throw new UnsupportedOperationException("Unknown database for match in protein MD5 ${proteinMd5}")
    }
}

def addSiteNodesDirect(locationSites, XMLStreamWriter writer) {
    writer.writeStartElement("sites")
    
    locationSites.each { siteMap ->
        writer.writeStartElement("site")
        description = siteMap.description ?: "null"
        writer.writeAttribute("description", description)
        writer.writeAttribute("numLocations", siteMap.numLocations.toString())
        
        writer.writeStartElement("site-locations")
        siteMap.siteLocations.each { siteLoc ->
            writer.writeStartElement("site-location")
            writer.writeAttribute("residue", siteLoc.residue)
            writer.writeAttribute("start", siteLoc.start.toString())
            writer.writeAttribute("end", siteLoc.end.toString())
            writer.writeEndElement()
        }
        writer.writeEndElement() // site-locations
        
        writer.writeEndElement() // site
    }
    
    writer.writeEndElement() // sites
}

def writeXrefDirect(seqData, XMLStreamWriter writer) {
    seqData.each { row ->
        writer.writeStartElement("xref")
        writer.writeAttribute("id", row.id)
        writer.writeAttribute("name", "${row.id} ${row.description}".trim())
        writer.writeEndElement()
    }
}

def fmtDefaultMatchNode(Map match) {
    return [
        evalue : match.evalue,
        score  : match.score
    ]
}

def fmtPantherMatchNode(Map match) {
    return [
        ac                 : match.treegrafter.subfamilyAccession,
        evalue             : match.evalue,
        "protein-class"    : match.treegrafter.proteinClass,
        "graft-point"      : match.treegrafter.graftPoint,
        "ancestral-node": match.treegrafter.ancestralNodeID,
        name               : match.signature.name,
        score              : match.score
    ]
}

def fmtPrintsMatchNode(Map match) {
    return [
        evalue    : match.evalue,
        graphscan : match.graphscan,
    ]
}

def fmtSuperfamilyMatchNode(Map match) {
    return [
        evalue : match.evalue
    ]
}

def fmtSignatureNode(Map match) {
    def signatureNodeAttributes = [ac: match.signature.accession]
    if (match.signature.name != null) {
        signatureNodeAttributes.name = match.signature.name
    }
    if (match.signature.desc != null) {
        signatureNodeAttributes.desc = match.signature.desc
    }
    if (match.signature.type != null) {
        signatureNodeAttributes.type = match.signature.type
    }
    return signatureNodeAttributes
}

// Keep all your existing location formatting methods:

def fmtDefaultLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        "hmm-start"    : loc.hmmStart,
        "hmm-end"      : loc.hmmEnd,
        "hmm-length"   : loc.hmmLength,
        "hmm-bounds"   : Location.getHmmBounds(loc.hmmBounds),
        evalue         : loc.evalue,
        score          : loc.score,
        "env-start"    : loc.envelopeStart,
        "env-end"      : loc.envelopeEnd
    ]
}

def fmtDefaultNoHbLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        "hmm-start"    : loc.hmmStart,
        "hmm-end"      : loc.hmmEnd,
        "hmm-length"   : loc.hmmLength,
        evalue         : loc.evalue,
        score          : loc.score,
        "env-start"    : loc.envelopeStart,
        "env-end"      : loc.envelopeEnd
    ]
}

def fmtMinimalistLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        score          : loc.score
    ]
}

def fmtCddLocationNode(Map match, Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        evalue         : match.evalue,
        score          : match.score,
    ]
}

def fmMinimalistLoctationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
    ]
}

def fmtMobidbLiteLocationNode(Map loc) {
    return [
        start              : loc.start,
        end                : loc.end,
        representative     : loc.representative,
        "sequence-feature" : loc.sequenceFeature,
    ]
}

def fmtPantherLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        "hmm-start"    : loc.hmmStart,
        "hmm-end"      : loc.hmmEnd,
        "hmm-length"   : loc.hmmLength,
        "hmm-bounds"   : Location.getHmmBounds(loc.hmmBounds),
        "env-start"    : loc.envelopeStart,
        "env-end"      : loc.envelopeEnd
    ]
}

def fmtPrintsLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        motifNumber    : loc.motifNumber,
        pvalue         : loc.pvalue,
        score          : loc.score
    ]
}

def fmtPrositePatternsLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        level          : loc.level,
    ]
}

def fmtSignalpLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        pvalue         : loc.score
    ]
}

def fmtSmartLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        evalue         : loc.evalue,
        score          : loc.score,
        "hmm-start"    : loc.hmmStart,
        "hmm-end"      : loc.hmmEnd,
        "hmm-length"   : loc.hmmLength,
        "hmm-bounds"   : Location.getHmmBounds(loc.hmmBounds)
    ]
}

def fmtSuperfamilyLocationNode(Map loc) {
    return [
        start          : loc.start,
        end            : loc.end,
        representative : loc.representative,
        evalue         : loc.evalue,
        "hmm-length"   : loc.hmmLength
    ]
}
