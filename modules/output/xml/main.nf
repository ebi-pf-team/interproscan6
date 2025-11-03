import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.core.*
import com.fasterxml.jackson.dataformat.xml.ser.ToXmlGenerator
import com.fasterxml.jackson.dataformat.xml.XmlFactory
import com.fasterxml.jackson.dataformat.xml.util.DefaultXmlPrettyPrinter
import java.io.StringWriter
import java.util.regex.Pattern
import javax.xml.namespace.QName
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.SeqDB

process WRITE_XML {
    label    'mem_low', 'time_long'
    executor 'local'

    input:
    val matches_files  // {query prot seq md5: {model acc: match}}
    val output_file
    val seq_db_file
    val nucleic
    val interproscan_version
    val db_releases

    exec:
    SeqDB db = new SeqDB(seq_db_file.toString())
    def xmlFactory = new XmlFactory()
    xmlFactory.configure(ToXmlGenerator.Feature.WRITE_XML_DECLARATION, true);
    def gen = (ToXmlGenerator) xmlFactory.createGenerator(new File(output_file), JsonEncoding.UTF8)
    gen.setPrettyPrinter(new DefaultXmlPrettyPrinter())
    gen.initGenerator()

    gen.setNextName(new QName("results"))
    gen.writeStartObject();
    gen.setNextIsAttribute(true);
    gen.writeFieldName("interproscan-version");
    gen.writeString(interproscan_version);
    gen.writeFieldName("interpro-version");
    gen.writeString(db_releases?.interpro?.version ?: "");
    Set<String> seenNucleicMd5s = new HashSet<>()

    matches_files.each { matchFile ->
        Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
        if (nucleic) {
            def (nucleicToProteinMd5, ntSeqDataMap, orfDataMap) =
                db.retrieveAllNucleicSequenceData(proteins.keySet() as List)

            nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                if (!seenNucleicMd5s.contains(nucleicMd5)) {
                    def seqData = ntSeqDataMap[nucleicMd5]
                    def protMd5ToOrfs = orfDataMap[nucleicMd5]
                    addNucleotideNode(gen, nucleicMd5, proteinMd5s, proteins, seqData, protMd5ToOrfs)
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
                addProteinNode(gen, proteinMd5, proteinMatches, seqData)
            }
        }
    }

    gen.writeEndObject()  // </results>
    gen.close()
}

def addNucleotideNode(ToXmlGenerator gen, String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, List ntSeqData, Map proteinMd5ToOrfs) {
    /* Write data for an input nucleic acid seq, and then the matches for its associated ORFs.
    <nucleotide-sequence>
        <sequence md5="" sequence </sequence>
        <xref id="id", name="id desc"/>
        <orf end="", start="", strand="">
            <protein>
                <sequence md5="" sequence </sequence>
                <xref id="id", name="id desc"/>
                <matches> <> <> </matches>
            </protein>
        </orf>
    */
    def SOURCE_NT_PATTERN = Pattern.compile(/^source=[^"]+\s+coords=(\d+)\.\.(\d+)\s+length=\d+\s+frame=(\d+)\s+desc=.*$/)

    gen.writeFieldName("nucleotide-sequence")
    gen.writeStartObject()

    gen.writeFieldName("sequence")
    gen.writeStartObject()
    gen.setNextIsAttribute(true);
    gen.writeStringField("md5", nucleicMd5)
    gen.setNextIsAttribute(false);
    gen.setNextIsUnwrapped(true);
    gen.writeStringField("", ntSeqData[0].sequence)
    gen.writeEndObject()

    writeXrefs(gen, ntSeqData)

    proteinMd5s.each { proteinMd5 ->
        // a proteinSeq MD5 may be associated with multiple nt seqs, only pull the data where the nt md5/seq is relevant
        def proteinSeqData = proteinMd5ToOrfs[proteinMd5]
        proteinSeqData.each { row ->
            def proteinSource = SOURCE_NT_PATTERN.matcher(row.description)
            assert proteinSource.matches()
            gen.writeFieldName("orf")
            gen.writeStartObject()
            gen.setNextIsAttribute(true);
            gen.writeNumberField("start", proteinSource.group(1) as int)
            gen.writeNumberField("start", proteinSource.group(2) as int)
            gen.writeStringField("start", proteinSource.group(3) as int < 4 ? "SENSE" : "ANTISENSE")  
            addProteinNode(gen, proteinMd5, proteinMatches[proteinMd5], proteinSeqData)
            gen.writeEndObject()
        }
    }

    gen.writeEndObject()
}

def addProteinNode(ToXmlGenerator gen, String proteinMd5, Map proteinMatches, List proteinSeqData) {  
    // Write data for a query protein sequence and its matches
    gen.writeFieldName("protein")
    gen.writeStartObject()

    // <sequence md5="...">...</sequence>
    gen.writeFieldName("sequence")
    gen.writeStartObject()
    gen.setNextIsAttribute(true);
    gen.writeStringField("md5", proteinMd5)
    gen.setNextIsAttribute(false);
    gen.setNextIsUnwrapped(true);
    gen.writeStringField("", proteinSeqData[0].sequence)
    gen.writeEndObject()

    // <xref id="..." name="..."/>
    writeXrefs(gen, proteinSeqData)

    // <matches>
    gen.writeFieldName("matches")
    gen.writeStartObject()
    proteinMatches.each { String modelAcc, Map match ->
        addMatchNode(gen, proteinMd5, match)
    }
    gen.writeEndObject()
    gen.writeEndObject()
}

def addMatchNode(ToXmlGenerator gen, String proteinMd5, Map match) {
    // Write an individual node representing a match. The structure is dependent on the application.
    String appl = (match.source == "InterPro-N"
            ? "InterPro-N" 
            : match.signature.signatureLibraryRelease.library
    ).toLowerCase()
    
    Map matchNodeAttributes = [source: match.source]
    switch (appl) {
        case "antifam":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-gene3d":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-funfam":
        case "funfam":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cdd":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "coils":
            break
        case "hamap":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "interpro-n":
            // No specific attributes for InterPro-N matches
            break
        case "mobidb lite":
        case "mobidb-lite":
        case "mobidb_lite":
            break
        case "ncbifam":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "panther":
            matchNodeAttributes = fmtPantherMatchNode(match)
            break
        case "pfam":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "phobius":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "pirsf":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "pirsr":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "prints":
            matchNodeAttributes = fmtPrintsMatchNode(match)
            break
        case "prosite patterns":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "prosite profiles":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "sfld":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "signalp":
            break
        case "smart":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "superfamily":
            matchNodeAttributes = fmtSuperfamilyMatchNode(match)
            break
        case "tmhmm":
        case "deeptmhmm":
        case "tmbed":
            break
        default:
            throw new UnsupportedOperationException("Unknown application '${appl}' for query protein with MD5 ${proteinMd5}")
    }

    // <match ...>
    gen.writeFieldName("match")
    gen.writeStartObject()
    gen.setNextIsAttribute(true)
    matchNodeAttributes.each { k, v ->
        if (v != null) {
            gen.writeStringField(k.toString(), v.toString())
        }
    }

    // <signature ...>
    gen.writeFieldName("signature")
    gen.writeStartObject()
    gen.setNextIsAttribute(true);
    def signatureNodeAttributes = fmtSignatureNode(match)
    signatureNodeAttributes.each { k, v ->
        if (v != null) {
            gen.writeStringField(k.toString(), v.toString())
        }
    }

    // <signature-library-release library="" version=""/>
    gen.writeFieldName("signature-library-release")
    gen.writeStartObject()
    gen.setNextIsAttribute(true);
    gen.writeStringField("library", match.signature.signatureLibraryRelease.library)
    gen.writeStringField("version", match.signature.signatureLibraryRelease.version)
    gen.writeEndObject()

    if (match.signature?.entry) {
        def entry = match.signature.entry
        addEntryNode(gen, entry)
    }

    gen.writeEndObject() // </signature>

    if (match.source != "InterPro-N" && appl == "panther") {
        gen.writeFieldName("model-ac")
        gen.writeString(match.treegrafter.subfamilyAccession ?: match.modelAccession)

        match.treegrafter.goXRefs.each { go ->
            gen.writeFieldName("go-xref")
            gen.writeStartObject()
            gen.setNextIsAttribute(true)
            gen.writeStringField("category", go.category)
            gen.writeStringField("db", go.databaseName)
            gen.writeStringField("id", go.id)
            gen.writeStringField("name", go.name)
            gen.writeEndObject()
        }
    } else {
        gen.writeFieldName("model-ac")
        gen.writeString(match.modelAccession)
    }

    addLocationNodes(gen, appl, proteinMd5, match)

    gen.writeEndObject() // </match>
}

def addEntryNode(ToXmlGenerator gen, Map entry) {
    /* Add info on the InterPro Entry the signature is integrated into. For example:
    <entry ac='IPR001584' desc='Integrase, catalytic core' name='Integrase_cat-core' type='Domain'>
        <go-xref category='BIOLOGICAL_PROCESS' db='GO' id='GO:0015074' name='DNA integration' />
        <pathway-xref db='MetaCyc' id='PWY-6955' name='lincomycin A biosynthesis' />
    </entry>
    */
    gen.writeFieldName("entry")
    gen.writeStartObject()
    gen.setNextIsAttribute(true)
    gen.writeStringField("ac", entry.accession)
    gen.writeStringField("desc", entry.description)
    gen.writeStringField("name", entry.name)
    gen.writeStringField("type", entry.type)

    if (entry.goXRefs) {
        entry.goXRefs.each { go ->
            gen.writeFieldName("go-xref")
            gen.writeStartObject()
            gen.setNextIsAttribute(true)
            gen.writeStringField("category", go.category)
            gen.writeStringField("db", go.databaseName)
            gen.writeStringField("id", go.id)
            gen.writeStringField("name", go.name)
            gen.writeEndObject()
        }
    }

    if (entry.pathwayXRefs) {
        entry.pathwayXRefs.each { pw ->
            gen.writeFieldName("pathway-xref")
            gen.writeStartObject()
            gen.setNextIsAttribute(true)
            gen.writeStringField("db", pw.databaseName)
            gen.writeStringField("id", pw.id)
            gen.writeStringField("name", pw.name)
            gen.writeEndObject()
        }
    }

    gen.writeEndObject()
}

// Formating the Match node

def fmtDefaultMatchNode(Map match) {
    return [
        evalue : match.evalue,
        score  : match.score,
        source : match.source
    ]
}

def fmtPantherMatchNode(Map match) {
    return [
        ac               : match.treegrafter.subfamilyAccession,
        evalue           : match.evalue,
        "protein-class"  : match.treegrafter.proteinClass,
        "graft-point"    : match.treegrafter.graftPoint,
        "ancestral-node" : match.treegrafter.ancestralNodeID,
        name             : match.signature.name,
        score            : match.score,
        source           : match.source
    ]
}

def fmtPrintsMatchNode(Map match) {
    return [
        evalue    : match.evalue,
        graphscan : match.graphscan,
        source    : match.source
    ]
}

def fmtSuperfamilyMatchNode(Map match) {
    return [
        evalue : match.evalue,
        source : match.source
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

// Formating and add Location nodes

def addLocationNodes(ToXmlGenerator gen, String appl, String proteinMd5, Map match) {
    gen.writeFieldName("locations")
    gen.writeStartObject()

    match.locations.each { loc ->
        def locationAttributes
        switch (appl) {
            case "antifam":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "cath-funfam":
            case "funfam":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "cath-gene3d":
            case "gene3d":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "cdd":
                locationAttributes = fmtCddLocationNode(match, loc)
                break
            case "coils":
                locationAttributes = fmMinimalistLocationNode(loc)
                break
            case "deeptmhmm":
                locationAttributes = []
                break
            case "hamap":
                locationAttributes = fmtMinimalistLocationWithScoreNode(loc)
                break
            case "interpro-n":
                locationAttributes = fmtMinimalistLocationWithScoreNode(loc)
                break
            case "mobidb lite":
            case "mobidb-lite":
            case "mobidb_lite":
                locationAttributes = fmtMobidbLiteLocationNode(loc)
                break
            case "ncbifam":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "panther":
                locationAttributes = fmtPantherLocationNode(loc)
                break
            case "pfam":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "phobius":
                locationAttributes = fmtMinimalistLocationWithScoreNode(loc)
                break
            case "pirsf":
                locationAttributes = fmtDefaultLocationNode(loc)
                break
            case "pirsr":
                locationAttributes = fmtDefaultNoHbLocationNode(loc)
                break
            case "prints":
                locationAttributes = fmtPrintsLocationNode(loc)
                break
            case "prosite patterns":
                locationAttributes = fmtPrositePatternsLocationNode(loc)
                break
            case "prosite profiles":
                locationAttributes = fmtMinimalistLocationWithScoreNode(loc)
                break
            case "sfld":
                locationAttributes = fmtDefaultNoHbLocationNode(loc)
                break
            case "signalp":
                locationAttributes = fmtSignalpLocationNode(loc)
                break
            case "smart":
                locationAttributes = fmtSmartLocationNode(loc)
                break
            case "superfamily":
                locationAttributes = fmtSuperfamilyLocationNode(loc)
                break
            case "tmhmm":
            case "deeptmhmm":
                locationAttributes = fmMinimalistLocationNode(loc)
                break
            case "tmbed":
                locationAttributes = fmMinimalistLocationNode(loc)
                break
            default:
                throw new UnsupportedOperationException("Unknown database for match ${matchId}")
        }

        gen.writeFieldName("location")
        gen.writeStartObject()
        gen.setNextIsAttribute(true)
        locationAttributes.each { k, v ->
            if (v != null) {
                if (v instanceof Number) {
                    gen.writeNumberField(k.toString(), ((Number) v).longValue())
                } else {
                    gen.writeStringField(k.toString(), v.toString())
                }
            }
        }

        if (loc.fragments && loc.fragments.size() > 0) {
            gen.writeFieldName("location-fragments")
            gen.writeStartObject()
            loc.fragments.each { frag ->
                gen.writeFieldName("fragment")
                gen.writeStartObject()
                gen.setNextIsAttribute(true)
                gen.writeNumberField("start", frag.start)
                gen.writeNumberField("end", frag.end)
                gen.writeStringField("dc-status", frag.dcStatus)
                gen.writeEndObject()
            }
            gen.writeEndObject()
        }

        if (appl in ["hamap", "prosite patterns", "prosite profiles"]) {
            gen.writeStringField("alignment", loc.targetAlignment ?: "")
            gen.writeStringField("cigar-alignment", loc.cigarAlignment ?: "")
        }

        if (loc.sites && loc.sites.size() > 0) {
            addSiteNodes(gen, loc.sites)
        }

        gen.writeEndObject()
    }

    gen.writeEndObject();
}

// formate location nodes

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

def fmtMinimalistLocationWithScoreNode(Map loc) {
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

def fmMinimalistLocationNode(Map loc) {
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

// add site nodes

def addSiteNodes(ToXmlGenerator gen, locationSites) {
    gen.writeFieldName("sites")
    gen.writeStartObject()
    locationSites.each { siteMap ->
        gen.writeFieldName("site")
        gen.writeStartObject()
        gen.setNextIsAttribute(true)
        gen.writeStringField("description", siteMap.description)
        gen.writeNumberField("numLocations", siteMap.numLocations)

        gen.writeFieldName("site-locations")
        gen.writeStartObject()
        siteMap.siteLocations.each { siteLoc ->
            gen.writeFieldName("site-location")
            gen.writeStartObject()
            gen.setNextIsAttribute(true)
            gen.writeStringField("residue", siteLoc.residue.toString())
            gen.writeNumberField("start", siteLoc.start)
            gen.writeNumberField("end", siteLoc.end)
            gen.writeEndObject()
        }
        gen.writeEndObject()
        gen.writeEndObject()
    }
    gen.writeEndObject()
}

def writeXrefs(ToXmlGenerator gen, List seqData) {
    seqData.each { row ->
        gen.writeFieldName("xref")
        gen.writeStartObject()
        gen.setNextIsAttribute(true);
        gen.writeStringField("id", row.id)
        gen.writeStringField("name", "${row.id} ${row.description ?: ''}".trim())
        gen.writeEndObject()
    }
}
