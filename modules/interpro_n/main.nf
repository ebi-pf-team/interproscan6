process PREPARE_INTERPRO_N {
    label    'mem_low'
    label    'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fasta)
    val batch_size    // max number of sequences in output fasta files
    val max_length
    val chunk_overlap

    output:
    tuple val(meta), path("sequences_*.tsv", arity: '1..*')

    exec:
    def sequences = uk.ac.ebi.interpro.FastaFile.chunkSequences(fasta, max_length, chunk_overlap)
    def fileIndex = 1
    def seqCount = 0
    def writer = task.workDir.resolve("sequences_${fileIndex}.tsv").newWriter("UTF-8")
    writer.writeLine("accession\tsequence")
    sequences.each { seq ->
        if (batch_size > 0 && seqCount >= batch_size) {
            writer.close()
            fileIndex += 1
            seqCount = 0
            writer = task.workDir.resolve("sequences_${fileIndex}.tsv").newWriter("UTF-8")
            writer.writeLine("accession\tsequence")
        }
        seq.chunks.eachWithIndex { chunk, idx ->
            writer.writeLine("${seq.id}_${idx + 1}\t${chunk}")
            seqCount += 1
        }
    }

    writer.close()
}

process RUN_INTERPRO_N_CPU {
    label     'mem_high'
    label     'time_medium'
    container 'interpro/interpro-n:1.0'

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
    cd /opt/interpro_n
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --output_dir \${WORKDIR}
    """
}

process RUN_INTERPRO_N_GPU {
    label     'mem_high'
    label     'time_short'
    label     'use_gpu'
    container 'interpro/interpro-n:1.0'

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
    cd /opt/interpro_n
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --output_dir \${WORKDIR}
    """
}

process PARSE_INTERPRO_N {
    label    'mem_low'
    label    'time_short'
    executor 'local'

    input:
    tuple val(meta), val(json_files), val(fasta)
    val applications
    val max_length
    val chunk_overlap
    val match_overlap

    output:
    tuple val(meta), path("interpro-n.json")

    exec:
    // Databases to report, given the user-requested applications
    def selectedApps = uk.ac.ebi.interpro.InterProN.selectDatabases(applications)

    // Sequence lengths, used to discard positions past the end of a chunk
    def lengths = uk.ac.ebi.interpro.FastaFile.parse(fasta).collectEntries { id, seq ->
        [(id): seq.length()]
    }

    // Parse all JSON files and group by sequence ID
    def jsonSlurper = new groovy.json.JsonSlurper()
    def resultsBySeq = [:].withDefault { [] }
    json_files.each { jsonFile ->
        def data = jsonSlurper.parse(jsonFile)
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
        chunks.sort { chunk -> chunk.id }

        def seqLength = lengths[seqId]
        assert seqLength != null

        // Collect all matched with global coordinates
        def allMatches = []
        chunks.each { chunk ->
            def offset = chunk.id * (max_length - chunk_overlap)
            // Last residue covered by this chunk
            def chunkEnd = Math.min(offset + max_length, seqLength)
            chunk.matches.each { match ->
                def sigLib = match.signature.signatureLibraryRelease
                def label = uk.ac.ebi.interpro.InterProN.toLabel(sigLib.library)
                if (label in selectedApps) {
                    sigLib.library = uk.ac.ebi.interpro.InterProN.SUPPORTED_DATABASES[label]
                } else {
                    return
                }
                assert match.locations.size() == 1

                def loc = match.locations[0]
                def start = loc.start + offset
                def end = Math.min(loc.end + offset, chunkEnd)
                if (start >= end) {
                    // Match on the end-of-sequence token only, or a single residue long
                    return
                }

                // Fragments beyond the last residue of the chunk only cover
                // the end-of-sequence token
                def fragments = loc["location-fragments"].findResults { frag ->
                    def fragStart = frag.start + offset
                    fragStart > chunkEnd
                        ? null
                        : [start: fragStart, end: Math.min(frag.end + offset, chunkEnd)]
                }
                if (fragments.any { frag -> frag.start >= frag.end }) {
                    // Single-residue fragments are spurious: discard the whole match
                    return
                }

                allMatches << [
                    signature: match.signature,
                    start: start,
                    end: end,
                    fragments: fragments,
                    score: loc.score as double  // stored as string in JSON
                ]
            }
        }

        // Group by signature
        def bySig = allMatches.groupBy { match ->
            def sig = match.signature
            def lib = sig.signatureLibraryRelease
            "${sig.accession}::${lib.library}::${lib.version}"
         }

        def merged = []
        bySig.each { _key, matches ->
            // Sort by descending score, then ascending start
            matches.sort { a, b ->
                def cmp = b.score <=> a.score
                (cmp != 0) ? cmp : (a.start <=> b.start)
            }

            def locations = []
            matches.each { match ->
                def matchLen = match.end - match.start + 1
                def mergedLoc = locations.find { other ->
                    def overlapStart = Math.max(other.start, match.start)
                    def overlapEnd = Math.min(other.end, match.end)
                    def overlapLength = overlapEnd - overlapStart + 1
                    if (overlapLength <= 0) return false

                    def otherLen = other.end - other.start + 1
                    overlapLength >= (match_overlap * otherLen) || overlapLength >= (match_overlap * matchLen)
                }

                if (mergedLoc) {
                    mergedLoc.start = Math.min(mergedLoc.start, match.start)
                    mergedLoc.end = Math.max(mergedLoc.end, match.end)
                    mergedLoc.fragments.addAll(match.fragments)
                    mergedLoc.fragments = mergeFragments(mergedLoc.fragments)
                } else {
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

            def library = new uk.ac.ebi.interpro.SignatureLibraryRelease(sigLib.library, sigLib.version)
            def signature = new uk.ac.ebi.interpro.Signature(accession, library)
            def match = new uk.ac.ebi.interpro.Match(accession, signature)
            // Override source set in Match constructor
            match.source = "InterPro-N"

            rawMatch.locations.each { loc ->
                // transform to LocationFragment objects
                def fragments = []
                if (loc.fragments.size() > 1) {
                    loc.fragments.sort { f -> f.start }
                    loc.fragments.eachWithIndex { f, i ->
                        def status
                        if (i == 0) {
                            status = "C_TERMINAL_DISC"
                        } else if (i + 1 < loc.fragments.size()) {
                            status = "NC_TERMINAL_DISC"
                        } else {
                            status = "N_TERMINAL_DISC"
                        }
                        fragments << new uk.ac.ebi.interpro.LocationFragment(f.start, f.end, status)
                    }
                } else {
                    def f = loc.fragments[0]
                    fragments << new uk.ac.ebi.interpro.LocationFragment(f.start, f.end, "CONTINUOUS")
                }

                match.addLocation(
                    new uk.ac.ebi.interpro.Location(loc.start, loc.end, loc.score, fragments)
                )
            }

            // use a namespaced key to avoid clashes with traditional applications
            def matchKey = "InterPro-N::${accession}".toString()
            results[seqId][matchKey] = match
        }
    }

    def filepath = task.workDir.resolve("interpro-n.json")
    filepath.text = groovy.json.JsonOutput.toJson(results)
}


def mergeFragments(fragments) {
    fragments.sort { f -> f.start }.inject([]) { merged, current ->
        if (merged.isEmpty()) {
            merged << current.clone()
            return merged
        }

        def last = merged[-1]
        if (current.start <= last.end + 1) {
            if (current.end > last.end) {
                last.end = current.end
            }
        } else {
            merged << current.clone()
        }
        merged
    }
}
