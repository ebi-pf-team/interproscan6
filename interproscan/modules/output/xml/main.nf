import com.fasterxml.jackson.databind.JsonNode
import com.fasterxml.jackson.databind.ObjectMapper
import com.fasterxml.jackson.databind.node.ObjectNode
import com.fasterxml.jackson.databind.SerializationFeature
import groovy.xml.MarkupBuilder
import java.time.format.DateTimeFormatter
import java.time.LocalDate
import java.io.StringWriter
import java.util.regex.Pattern

process WRITE_XML_OUTPUT {
    label 'local', 'ips6_container'

    input:
    val matchesFiles  // {query prot seq md5: {model acc: match}}
    val outputPath
    val seqDbPath
    val nucleic
    val ips6Version

    exec:
    def writer = new StringWriter()
    def xml = new MarkupBuilder(writer)
    // set the correct encoding so symbols are formatted correctly in the final output
    xml.setEscapeAttributes(false) // Prevent escaping attributes
    xml.setEscapeText(false)       // Prevent escaping text
    def analysisType = nucleic ? "nucleotide-sequence-matches" : "protein-matches"
    def jacksonMapper = new ObjectMapper().enable(SerializationFeature.INDENT_OUTPUT)
    def seqDatabase = new SeqDatabase(seqDbPath.toString())

    xml."$analysisType"("interproscan-version": ips6Version) {
        matchesFiles.each { matchFile ->
            if (nucleic) {
                proteinMatches = JsonReader.load(matchFile.toString(), jacksonMapper)
                nucleicToProteinMd5 = seqDatabase.groupProteins(proteinMatches)
                nucleicToProteinMd5.each { String nucleicMd5, Set<String> proteinMd5s ->
                    addNucleotideNode(nucleicMd5, proteinMd5s, proteinMatches, xml)
                }
            } else {
                JsonReader.streamJson(matchFile.toString(), jacksonMapper) { String proteinMd5, JsonNode matchesNode ->
                    addProteinNodes(proteinMd5, matchesNode, xml, seqDatabase)
                }
            }
        }
    }

    def outputFilePath = "${outputPath}.ips6.xml"
    new File(outputFilePath).text = writer.toString()
}

def addNucleotideNode(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, def xml, SeqDatabase seqDatabase) {
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

    // 1. <nt-seq> <sequence md5="" "<seq>" </sequence>
    ntSeqData = seqDatabase.getSeqData(nucleicMd5, true)
    String sequence = ntSeqData[0].split('\t')[-1]
    xml."nucleotideNode" {
        sequence(md5: nucleicMd5, sequence)

        // 2. <xref id="id" name="id desc"/>
        writeXref(seqData, xml)

        // 3. <orf end="", start="", strand="">
        proteinMd5s.each { proteinMd5 ->
            // a proteinSeq/Md5 may be associated with multiple nt md5s/seq, only pull the data where the nt md5/seq is relevant
            proteinSeqData = seqDatabase.getSeqData(proteinMd5, false, nucleicMd5)
            proteinSeqData.each { row ->
                def proteinDesc = row.split("\t")[1]
                def proteinSource = SOURCE_NT_PATTERN.matcher(proteinDesc)
                assert proteinSource.matches()
                orf(
                    start  : proteinSource.group(1) as int,
                    end    : proteinSource.group(2) as int,
                    strand : proteinSource.group(3) as int < 4 ? "SENSE" : "ANTISENSE"
                ) {
                    // 4. <protein> ... <\protein>
                    addProteinNodes(proteinMd5, proteinMatches[proteinMd5], xml, seqDatabase)
                }
            }
        }
    }
}

def addProteinNodes (String proteinMd5, JsonNode matchesNode, def xml, SeqDatabase seqDatabase) {
    /* Write data for a query protein sequence and its matches:
    <protein>
        <sequence md5="" sequence </sequence>
        <xref id="id", name="id desc"/>
        <matches> <> <> </matches>
    </protein>
    There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */
    xml.protein {
        // 1. <sequence md5="" sequence </sequence>
        proteinSeqData = seqDatabase.getSeqData(proteinMd5, false)
        String sequence = proteinSeqData[0].split('\t')[-1]
        xml.sequence(md5: proteinMd5, sequence)

        // 2. <xref id="id", name="id desc"/>
        writeXref(proteinSeqData, xml)

        // 3. <matches> <m> <m> <m> </matches>
        matches {
            matchesNode.fields().each { Map.Entry<String, JsonNode> entry ->
                JsonNode match = entry.value
                addMatchNode(proteinMd5, match, xml)
            }
        }
    }
}

def addMatchNode(String proteinMd5, JsonNode match, def xml) {
    // Write an individual node representing a match. The structure is dependent on the memberDB.
    String memberDB = JsonReader.asString(match.get("signature").get("signatureLibraryRelease").get("library")).toLowerCase() ?: ""
    switch (memberDB) {
        case "antifam":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-gene3d":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-funfam":
        case "funfam":  // use groovy case fall to allow multiple options
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cdd":
            matchNodeName = "CDD"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "coils":
            matchNodeName = "COILS"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "hamap":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "mobidb-lite":
        case "mobidb_lite":  // use groovy case fall to allow multiple options
            matchNodeName = "Mobidb-Lite"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "ncbifam":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "panther":
            matchNodeName = "PANTHER"
            matchNodeAttributes = fmtPantherMatchNode(match)
            break
        case "pfam":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "phobius":
            matchNodeName = "Phobius"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "pirsf":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "pirsr":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "prints":
            matchNodeName = "PRINTS"
            matchNodeAttributes = fmtPrintsMatchNode(match)
            break
        case "prosite patterns":
            matchNodeName = "PROSITE-patterns"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "prosite profiles":
            matchNodeName = "PROSITE-profiles"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "sfld":
            matchNodeName = "SFLD"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "signalp":
            matchNodeName = "SignalP"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "smart":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "superfamily":
            matchNodeName = "hmmer3"
            matchNodeAttributes = fmtSuperfamilyMatchNode(match)
            break
        default:
            throw new UnsupportedOperationException("Unknown database '${memberDB}' for query protein with MD5 ${proteinMd5}")
    }

    def signatureNodeAttributes = fmtSignatureNode(match)

    xml."$matchNodeName-match"(matchNodeAttributes) {
        xml.signature(signatureNodeAttributes) {
            xml.signatureLibraryRelease {
                xml.library(match.get("signature").get("signatureLibraryRelease").get("library"))
                xml.version(match.get("signature").get("signatureLibraryRelease").get("version"))
            }
            // GO TERMS AND PATHWAYS
        }
        xml."model-ac"(memberDB == "panther" ? match.treegrafter.subfamilyAccession : match.modelAccession)

        addLocationNodes(matchNodeName, memberDB, proteinMd5, match, xml)
    }
}

// Formating the Match node

def fmtDefaultMatchNode(JsonNode source) {
    return [
        evalue : source.get("evalue"),
        score  : source.get("score")
    ]
}

def fmtPantherMatchNode(JsonNode source) {
    return [
        ac            : JsonReader.asString(source.get("treegrafter").get("subfamilyAccession")),
        evalue        : source.get("evalue"),
        "graft-point" : JsonReader.asString(source.get("treegrafter").get("graftPoint")),
        name          : JsonReader.asString(source.get("signature").get("name")),
        score         : source.get("score")
    ]
}

def fmtPrintsMatchNode(JsonNode source) {
    return [
        evalue    : source.get("evalue"),
        graphscan : JsonReader.asString(source.get("graphScan")),
    ]
}

def fmtSuperfamilyMatchNode(JsonNode source) {
    return [
        evalue : source.get("evalue"),
    ]
}

// Formating the Signature node

def fmtSignatureNode(JsonNode match) {
    def name = JsonReader.asString(match.get("signature").get("name"))
    def desc = JsonReader.asString(match.get("signature").get("description"))
    return [ac: JsonReader.asString(match.get("signature").get("accession"))].with {
       if (name) name = name
       if (desc) desc = desc
       it
   }
}

// Formating and add Location nodes

def addLocationNodes(String matchNodeName, String memberDB, String proteinMd5, JsonNode match, def xml) {
    xml.locations {
        match.get("locations").each { loc ->
            def locationAttributes
            switch (memberDB) {
                case "antifam":
                    locationAttributes = fmtDefaultLocationNode(loc)
                    break
                case "cath-funfam":
                case "cath-gene3d":
                    locationAttributes = fmtDefaultLocationNode(loc)
                    break
                case "cdd":
                    locationAttributes = fmtCddLocationNode(loc, match)
                    break
                case "coils":
                    locationAttributes = []
                    break
                case "deeptmhmm":
                    locationAttributes = []
                    break
                case "funfam":
                    locationAttributes = fmtDefaultLocationNode(loc)
                    break
                case "hamap":
                    locationAttributes = fmtMinimalistLocationNode(loc)
                    break
                case "mobidb-lite":
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
                    locationAttributes = []
                    break
                case "pirsf":
                case "pirsr":
                    locationAttributes = fmtDefaultLocationNode(loc)
                    break
                case "prints":
                    locationAttributes = fmtPrintsLocationNode(loc)
                    break
                case "prosite patterns":
                    locationAttributes = fmtPrositePatternsLocationNode(loc)
                    break
                case "prosite profiles":
                    locationAttributes = fmtMinimalistLocationNode(loc)
                    break
                case "sfld":
                    locationAttributes = fmtDefaultLocationNode(loc)
                    break
                case "signalp":
                    locationAttributes = fmtMinimalistLocationNode(loc)
                    break
                case "smart":
                    locationAttributes = fmtSmartLocationNode(loc)
                    break
                case "superfamily":
                    locationAttributes = fmtSuperfamilyLocationNode(loc)
                    break
                default:
                    throw new UnsupportedOperationException("Unknown database for match ${matchId}")
            }

            xml.location(locationAttributes) {
                if (loc.has("location-fragments") && loc.get("location-fragments").size() > 0) {
                    xml."$matchNodeName-fragment" {
                        loc.get("location-fragments").each { frag ->
                            start(frag.get("start"))
                            end(frag.get("end"))
                            dcStatus(frag.get("dcStatus"))
                        }
                    }
                }

                if (memberDB in ["hamap", "prosite patterns", "prosite profiles"]) {
                    xml.alignment(JsonReader.asString(loc.get("targetAlignment") ?: ""))
                }
                if (loc.has("sites")) {
                    if (loc.get('sites').size() > 0) {
                        xml.addSiteNodes(loc.get("sites"), memberDB, xml)
                    }
                }
            }
        }
    }
}

// formate location nodes

def fmtDefaultLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        hmmStart       : loc.get("hmmStart"),
        hmmEnd         : loc.get("hmmEnd"),
        hmmLength      : loc.get("hmmLength"),
        hmmBounds      : loc.get("hmmBounds"),
        evalue         : loc.get("evalue"),
        score          : loc.get("score"),
        envelopeStart  : loc.get("envelopeStart"),
        envelopeEnd    : loc.get("envelopeEnd")
    ]
}

def fmtMinimalistLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        score          : loc.get("score"),
    ]
}

def fmtCddLocationNode(match, loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        evalue         : match.get("evalue"),
        score          : match.get("score"),
    ]
}

def fmtMobidbLiteLocationNode(Map loc) {
    return [
        start              : loc.get("start"),
        end                : loc.get("end"),
        representative     : loc.get("representative"),
        "sequence-feature" : loc.get("sequenceFeature"),
    ]
}

def fmtPantherLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        hmmStart       : loc.get("hmmStart"),
        hmmEnd         : loc.get("hmmEnd"),
        hmmLength      : loc.get("hmmLength"),
        hmmBounds      : loc.get("hmmBounds"),
        envelopeStart  : loc.get("envelopeStart"),
        envelopeEnd    : loc.get("envelopeEnd")
    ]
}

def fmtPrintsLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        motifNumber    : loc.get("motifNumber"),
        pvalue         : loc.get("pvalue"),
        score          : loc.get("score")
    ]
}

def fmtPrositePatternsLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        level          : loc.get("level"),
    ]
}

def fmtSmartLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        evalue         : loc.get("evalue"),
        score          : loc.get("score"),
        hmmStart       : loc.get("hmmStart"),
        hmmEnd         : loc.get("hmmEnd"),
        hmmLength      : loc.get("hmmLength"),
        hmmBounds      : loc.get("hmmBounds")
    ]
}

def fmtSuperfamilyLocationNode(Map loc) {
    return [
        start          : loc.get("start"),
        end            : loc.get("end"),
        representative : loc.get("representative"),
        evalue         : loc.get("evalue"),
        hmmLength      : loc.get("hmmLength")
    ]
}

// add site nodes

// TODO: change to json node handling
def addSiteNodes(locationSites, memberDB, xml) {
    xml."sites" {
        locationSites.each { siteObj ->
            xml."$matchNodeName-site"(description: JsonReader.asString(siteObj.get("description")), numLocations: JsonReader.getString(siteObj.get("numLocations"))) {
                if(siteObj.group){ group(siteObj.get("group")) }
                if(siteObj.label){ label(siteObj.get("label")) }
                if (memberDB != "cdd") {
                    hmmStart(siteObj.get("hmmStart"))
                    hmmEnd(siteObj.get("hmmEnd"))
                }
                "site-locations" {
                    siteObj.siteLocations.each { siteLoc ->
                        xml."site-location"(
                            residue : siteLoc.get("residue"),
                            start   : siteLoc.get("start"),
                            end     : siteLoc.get("end")
                        )
                    }
                }
            }
        }
    }
}

def writeXref(seqData, xml) {
    // <xref id="id" name="id desc"/>
    seqData.each { row ->
        row = row.split('\t')
        seqId = row[0]
        seqDesc = row[1]
        xml.xref(id: seqId, name: "$seqId $seqDesc")
    }
}
