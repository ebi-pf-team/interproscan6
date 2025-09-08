import com.fasterxml.jackson.databind.ObjectMapper
import groovy.xml.StreamingMarkupBuilder
import java.io.FileWriter
import java.util.regex.Pattern

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
    def fileWriter = new FileWriter(output_file)
    def BATCH_SIZE = 10000

    try {
        Set<String> seenNucleicMd5s = new HashSet<>()
        def builder = new StreamingMarkupBuilder()
        builder.encoding = 'UTF-8'

        def xmlContent = builder.bind { // mkp = Markup Builder; delegate = the XML builder context
            mkp.xmlDeclaration()
            results("interproscan-version": interproscan_version, "interpro-version": db_releases?.interpro?.version) {
                matches_files.each { matchFile ->
                    Map proteins = new ObjectMapper().readValue(new File(matchFile.toString()), Map)
                    if (nucleic) {
                        processNucleotidesBulk(db, proteins, seenNucleicMd5s, delegate, processedCount, BATCH_SIZE)
                    } else {
                        def proteinMd5List = proteins.keySet().toList()
                        Map<String, List> seqData = db.proteinMd5sToProteinSeqs(proteinMd5List)

                        for (Map.Entry<String, Map> entry : proteins.entrySet()) {
                            String proteinMd5 = entry.key
                            Map proteinMatches = entry.value
                            addProteinNodes(proteinMd5, proteinMatches, seqData[proteinMd5], delegate)
                        }
                    }
                }
            }
        }

        fileWriter << xmlContent
    } finally {
        fileWriter.close()
        db.close()
    }
}

def processNucleotidesBulk(SeqDBQuery db, Map proteins, Set seenNucleicMd5s, def xml, int processedCount, int batchSize) {
    Set<String> allProteinMd5s = proteins.keySet().toSet()
    Map<String, Set<String>> nucleicToProteinMd5 = db.groupProteinsBulk(allProteinMd5s)
    
    // Filter to only unseen nucleotide MD5s
    List<String> newNucleicMd5s = nucleicToProteinMd5.keySet().findAll { !seenNucleicMd5s.contains(it) }
    seenNucleicMd5s.addAll(newNucleicMd5s)
    if (newNucleicMd5s.isEmpty()) {
        return processedCount
    }

    // Get all nucleotide sequences
    Map<String, List> nucleotideSeqData = db.nucleicMd5sToNucleicSeqs(newNucleicMd5s)

    // Get all ORF sequences for relevant protein/nucleotide combinations
    List<String> relevantProteinMd5s = []
    newNucleicMd5s.each { nucleicMd5 ->
        relevantProteinMd5s.addAll(nucleicToProteinMd5[nucleicMd5])
    }
    Map<String, Map<String, List>> orfSeqData = db.getOrfSeqsBulk(relevantProteinMd5s, newNucleicMd5s)

    // Process the results
    newNucleicMd5s.each { String nucleicMd5 ->
        def ntSeqData = nucleotideSeqData[nucleicMd5]
        if (!ntSeqData) return

        Set<String> proteinMd5sForNucleic = nucleicToProteinMd5[nucleicMd5]

        addNucleotideNodeBulk(nucleicMd5, proteinMd5sForNucleic, proteins, ntSeqData, orfSeqData, xml)
    }
    
    return processedCount
}

def addNucleotideNode(String nucleicMd5, Set<String> proteinMd5s, Map proteinMatches, List ntSeqData, Map orfSeqData,def xml) {
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
    String sequence = ntSeqData[0].sequence
    xml."nucleotideNode" {
        xml.sequence(md5: nucleicMd5, sequence)

        // 2. <xref id="id" name="id desc"/>
        writeXref(ntSeqData, xml)

        // 3. <orf end="", start="", strand="">
        proteinMd5s.each { proteinMd5 ->
            // a proteinSeq MD5 may be associated with multiple nt seqs, only pull the data where the nt md5/seq is relevant
            proteinSeqData = orfSeqData[proteinMd5]?.get(nucleicMd5)
            if (proteinSeqData) {
                proteinSeqData.each { row ->
                    def proteinSource = SOURCE_NT_PATTERN.matcher(row.description)
                    assert proteinSource.matches()
                    xml.orf([
                        start  : proteinSource.group(1) as int,
                        end    : proteinSource.group(2) as int,
                        strand : proteinSource.group(3) as int < 4 ? "SENSE" : "ANTISENSE"
                    ]) {
                        // 4. <protein> ... <\protein>
                        addProteinNodes(proteinMd5, proteinMatches[proteinMd5], proteinSeqData[proteinMd5], xml)
                    }
                }
            }
        }
    }
}

def addProteinNodes (String proteinMd5, Map proteinMatches, List proteinSeqData, def xml) {
    /* Write data for a query protein sequence and its matches:
    <protein>
        <sequence md5="" sequence </sequence>
        <xref id="id", name="id desc"/>
        <matches> <> <> </matches>
    </protein>
    There may be multiple seqIds and desc for the same sequence/md5, use the first entry to get the seq. */
    String sequence = proteinSeqData[0].sequence
    xml.protein {
        // 1. <sequence md5="" sequence </sequence>
        xml.sequence(md5: proteinMd5, sequence)

        // 2. <xref id="id", name="id desc"/>
        writeXref(proteinSeqData, xml)

        // 3. <matches> <m> <m> <m> </matches>
        matches {
            proteinMatches.each { String modelAcc, Map match ->
                addMatchNode(proteinMd5, match, xml)
            }
        }
    }
}

def addMatchNode(String proteinMd5, Map match, def xml) {
    // Write an individual node representing a match. The structure is dependent on the memberDB.
    String memberDB = match.signature.signatureLibraryRelease.library.toLowerCase() ?: ""

    // Define the name for the match node and it's attributes
    switch (memberDB) {
        case "antifam":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-gene3d":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cath-funfam":
        case "funfam":  // use groovy case fall to allow multiple options
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "cdd":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "coils":
            matchNodeAttributes = null
            break
        case "hamap":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "mobidb lite":
        case "mobidb-lite":
        case "mobidb_lite":  // use groovy case fall to allow multiple options
            matchNodeAttributes = null
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
            matchNodeAttributes =  null
            break
        case "smart":
            matchNodeAttributes = fmtDefaultMatchNode(match)
            break
        case "superfamily":
            matchNodeAttributes = fmtSuperfamilyMatchNode(match)
            break
        case "tmhmm":
        case "deeptmhmm":
            matchNodeAttributes = null
            break
        case "tmbed":
            matchNodeAttributes = null
            break
        default:
            throw new UnsupportedOperationException("Unknown database '${memberDB}' for query protein with MD5 ${proteinMd5}")
    }

    def signatureNodeAttributes = fmtSignatureNode(match)

    xml."match"(matchNodeAttributes) {
        xml.signature(signatureNodeAttributes) {
            if (match.signature.entry && match.signature.entry != null) {
                addEntryNode(match.signature.entry, xml)
            }
            xml."signature-library-release"(
                library: match.signature.signatureLibraryRelease.library,
                version: match.signature.signatureLibraryRelease.version
            )
        }

        if (memberDB == "panther") {
            xml."model-ac"(match.treegrafter.subfamilyAccession ?: match.modelAccession)
            match.treegrafter.goXRefs.each { goXref ->
                xml."go-xref"(
                    category: goXref.category,
                    db: goXref.databaseName,
                    id: goXref.id,
                    name: goXref.name
                )
            }
        } else {
            xml."model-ac"(match.modelAccession)
        }

        addLocationNodes(memberDB, proteinMd5, match, xml)
    }
}

def addEntryNode(Map entry, def xml) {
    /* Add info on the InterPro Entry the signature is integrated into. For example:
    <entry ac='IPR001584' desc='Integrase, catalytic core' name='Integrase_cat-core' type='Domain'>
        <go-xref category='BIOLOGICAL_PROCESS' db='GO' id='GO:0015074' name='DNA integration' />
        <pathway-xref db='MetaCyc' id='PWY-6955' name='lincomycin A biosynthesis' />
    </entry>
    */
    xml.entry(
        ac: entry.accession,
        desc: entry.description,
        name: entry.name,
        type: entry.type
    ) {
        if (entry.goXRefs != null) {
            entry.goXRefs.each { goXref ->
                xml."go-xref"(
                    category: goXref.category,
                    db: goXref.databaseName,
                    id: goXref.id,
                    name: goXref.name
                )
            }
        }
        if (entry.pathwayXRefs != null) {
            entry.pathwayXRefs.each { pathwayXref ->
                xml."pathway-xref"(
                    db: pathwayXref.databaseName,
                    id: pathwayXref.id,
                    name: pathwayXref.name
                )
            }
        }
    }
}

// Formating the Match node

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

// Formating and add Location nodes

def addLocationNodes(String memberDB, String proteinMd5, Map match, def xml) {
    xml.locations {
        match.locations.each { loc ->
            def locationAttributes
            switch (memberDB) {
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
                    locationAttributes = fmMinimalistLoctationNode(loc)
                    break
                case "deeptmhmm":
                    locationAttributes = []
                    break
                case "hamap":
                    locationAttributes = fmtMinimalistLocationNode(loc)
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
                    locationAttributes = fmtMinimalistLocationNode(loc)
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
                    locationAttributes = fmtMinimalistLocationNode(loc)
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
                    locationAttributes = fmMinimalistLoctationNode(loc)
                    break
                case "tmbed":
                    locationAttributes = fmMinimalistLoctationNode(loc)
                    break
                default:
                    throw new UnsupportedOperationException("Unknown database for match ${matchId}")
            }

            xml.location(locationAttributes) {
                if (loc.containsKey("fragments") && loc["fragments"].size() > 0) {
                    xml."location-fragments" {
                        loc.fragments.each { frag ->
                            xml.fragment([
                                start      : frag.start,
                                end        : frag.end,
                                "dc-status": frag.dcStatus
                            ])
                        }
                    }
                }
                if (memberDB in ["hamap", "prosite patterns", "prosite profiles"]) {
                    xml.alignment(loc.targetAlignment ?: "")
                    xml."cigar-alignment"(loc.cigarAlignment ?: "")
                }
                if (loc.containsKey("sites") && loc.sites.size() > 0) {
                    addSiteNodes(loc.sites, memberDB, xml)
                }
            }
        }
    }
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

// add site nodes

def addSiteNodes(locationSites, memberDB, xml) {
    xml."sites" {
        locationSites.each { siteMap ->
            xml."site"(description: siteMap.description, numLocations: siteMap.numLocations) {
                "site-locations" {
                    siteMap.siteLocations.each { siteLoc ->
                        xml."site-location"(
                            residue : siteLoc.residue,
                            start   : siteLoc.start,
                            end     : siteLoc.end
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
        xml.xref(id: row.id, name: "${row.id} ${row.description}".trim())
    }
}
