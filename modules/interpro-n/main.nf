// codenarc-disable AllowedDirectivesRule
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.LocationFragment
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease

def MAX_LENGTH    = 2047  // Maximum sequence length
def CHUNK_OVERLAP = 1000  // Number of overlapping residues between consecutive chunks
def MATCH_OVERLAP = 0.25  // Fractional overlap threshold for merging matches
def DATABASES = [         // Databases supported by InterPro-N
    cathgene3d     : "CATH-Gene3D",
    cdd            : "CDD",
    hamap          : "HAMAP",
    ncbifam        : "NCBIFAM",
    panther        : "PANTHER",
    pfam           : "Pfam",
    pirsf          : "PIRSF",
    prints         : "PRINTS",
    prositeprofiles: "PROSITE profiles",
    prositepatterns: "PROSITE patterns",
    sfld           : "SFLD",
    smart          : "SMART",
    ssf            : "SUPERFAMILY",
]

process PREPARE_INTERPRO_N {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fasta)
    val batch_size    // max number of sequences in output fasta files

    output:
    tuple val(meta), path("sequences_*.tsv", arity: '1..*')

    exec:
    def sequences = FastaFile.chunkSequences(fasta.toString(), MAX_LENGTH, CHUNK_OVERLAP)
    def fileIndex = 1
    def seqCount = 0
    def writer = null
    def openNewFile = {
        if (writer) writer.close()
        def filePath = task.workDir.resolve("sequences_${fileIndex}.tsv").toFile()
        writer = filePath.newWriter("UTF-8")
        writer.writeLine("accession\tsequence")
        fileIndex++
        seqCount = 0
        return writer
    }

    writer = openNewFile()
    sequences.each { seq ->
        if (batch_size > 0 && seqCount >= batch_size) {
            writer = openNewFile()
        }
        seq.chunks.eachWithIndex { chunk, idx ->
            writer.writeLine("${seq.id}_${idx + 1}\t${chunk}")
            seqCount++
        }
    }

    if (writer) writer.close()
}

process RUN_INTERPRO_N_CPU {
    label 'mem_high', 'time_medium', 'interpro_n_container'

    input:
    tuple val(meta), path(tsv)
    path ck_dir
    val batch_size

    output:
    tuple val(meta), path("*json", arity: '1..*')

    script:
    """
    WORKDIR=\$(pwd)
    INPUT=\$(readlink -e ${tsv})
    CHECKPOINT=\$(readlink -e ${ck_dir})
    cd /workdir
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --output_dir \${WORKDIR}
    """
}

process RUN_INTERPRO_N_GPU {
    label 'mem_high', 'time_short', 'interpro_n_container', 'use_gpu'

    input:
    tuple val(meta), path(tsv)
    path ck_dir
    val batch_size

    output:
    tuple val(meta), path("*json", arity: '1..*')

    script:
    """
    WORKDIR=\$(pwd)
    INPUT=\$(readlink -e ${tsv})
    CHECKPOINT=\$(readlink -e ${ck_dir})
    cd /workdir
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --output_dir \${WORKDIR}
    """
}

process PARSE_INTERPRO_N {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(json_files)
    val applications

    output:
    tuple val(meta), path("interpro-n.json")

    exec:
    // User-requested applications
    def requested = applications.collect { name ->
        switch(name) {
            case "superfamily": return "ssf"
            default           : return name
        }
    } as Set

    // Compute overlap
    def overlappingApps = DATABASES.keySet().intersect(requested)

    // If none overlap, fallback to all supported applications
    def selectedApps = overlappingApps ?: DATABASES.keySet()

    // Parse all JSON files and group by sequence ID
    def jsonSlurper = new JsonSlurper()
    def resultsBySeq = [:].withDefault { [] }
    json_files.each { jsonFile ->
        def file = new File(jsonFile.toString())
        def data = jsonSlurper.parse(file)
        data.each { item ->
            def m = item.accession =~ /^(.+)_([0-9]+)$/
            assert m.size() == 1
            def seqId = m[0][1]
            def chunkId = m[0][2].toInteger()
            resultsBySeq[seqId] << [id: chunkId - 1, matches: item.match]
        }
    }

    // Merge matches for each full sequence
    def results = [:]
    resultsBySeq.each { seqId, chunks ->
        chunks.sort { it.id }

        // Collect all matched with global coordinates
        def allMatches = []
        chunks.each { chunk ->
            def offset = chunk.id * (MAX_LENGTH - CHUNK_OVERLAP)
            chunk.matches.each { match ->
                def sigLib = match.signature.signatureLibraryRelease
                def sigLibName = sigLib.library.toLowerCase().replaceAll(/[-\s]/, "")
                if (sigLibName in selectedApps) {
                    sigLib.library = DATABASES[sigLibName]
                } else {
                    return
                }
                assert match.locations.size() == 1

                def loc = match.locations[0]   
                allMatches << [
                    signature: match.signature,
                    start: loc.start + offset,
                    end: loc.end + offset,
                    fragments: loc["location-fragments"].collect { frag ->
                        [start: frag.start + offset, end: frag.end + offset]
                    },
                    score: loc.score as double  // stored as string in JSON
                ]
            }
        }

        // Group by signature
        def bySig = allMatches.groupBy {
            def sig = it.signature
            def lib = sig.signatureLibraryRelease
            "${sig.accession}::${lib.library}::${lib.version}"
         }

        def merged = []
        bySig.each { sigKey, matches ->
            // Sort by descending score, then ascending start
            matches.sort { a, b ->
                int cmp = b.score <=> a.score
                (cmp != 0) ? cmp : (a.start <=> b.start)
            }

            def locations = []
            matches.each { match ->
                boolean isMerged = false

                for (other in locations) {
                    def overlapStart = Math.max(other.start, match.start)
                    def overlapEnd = Math.min(other.end, match.end)
                    def overlapLength = overlapEnd - overlapStart + 1
                    if (overlapLength <= 0) continue

                    def lenA = other.end - other.start + 1
                    def lenB = match.end - match.start + 1
                    def fracA = overlapLength / lenA
                    def fracB = overlapLength / lenB

                    if (fracA >= MATCH_OVERLAP || fracB >= MATCH_OVERLAP) {
                        other.start = Math.min(other.start, match.start)
                        other.end = Math.max(other.end, match.end)
                        other.fragments.addAll(match.fragments)
                        other.fragments = mergeFragments(other.fragments)
                        isMerged = true
                        break
                    }
                }
                if (!isMerged) {
                    match.fragments = mergeFragments(match.fragments)
                    locations << match
                }
            }

            merged << [
                signature: locations[0].signature,
                locations: locations
            ]
        }

        merged.each { rawMatch ->           
            results.computeIfAbsent(seqId) { [:] }
            
            def sig = rawMatch.signature
            def sigLib = sig.signatureLibraryRelease
            def accession = sig.accession

            SignatureLibraryRelease library = new SignatureLibraryRelease(sigLib.library, sigLib.version)
            Signature signature = new Signature(accession, library)
            Match match = new Match(accession, signature)
            // Override source set in Match constructor
            match.source = "InterPro-N"

            rawMatch.locations.each { loc ->
                // transform to LocationFragment objects
                List<LocationFragment> fragments = []
                if (loc.fragments.size() > 1) {
                    loc.fragments.sort { it.start }
                    loc.fragments.eachWithIndex { f, i ->
                        String status
                        if (i == 0) {
                            status = "C_TERMINAL_DISC"
                        } else if (i + 1 < loc.fragments.size()) {
                            status = "NC_TERMINAL_DISC"
                        } else {
                            status = "N_TERMINAL_DISC"
                        }
                        fragments << new LocationFragment(f.start, f.end, status)
                    }
                } else {
                    def f = loc.fragments[0]
                    fragments << new LocationFragment(f.start, f.end, "CONTINUOUS")
                }

                match.addLocation(
                    new Location(loc.start, loc.end, loc.score, fragments)
                )
            }

            // use a namespaced key to avoid clashes with traditional applications
            def matchKey = "InterPro-N::${accession}".toString()
            results[seqId][matchKey] = match
        }
    }

    def json = JsonOutput.toJson(results)
    def outputFilePath = task.workDir.resolve("interpro-n.json")
    new File(outputFilePath.toString()).write(json)
}


def mergeFragments(fragments) {
    def sorted = fragments.sort { it.start }
    def merged = [sorted[0].clone()]
    for (int i = 1; i < sorted.size(); i++) {
        def current = sorted[i]
        def last = merged[-1]
        if (current.start <= last.end + 1) {
            last.end = Math.max(last.end, current.end)
        } else {
            merged << current.clone()
        }
    }
    return merged
}