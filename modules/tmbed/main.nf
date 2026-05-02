process PREPARE_TMBED {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fasta)
    val chunk_size    // max number of aa per chunked sequence
    val chunk_overlap // overlap between two chunks of a sequence
    val batch_size    // max number of sequences in output fasta files

    output:
    tuple val(meta), path("tmbed_*.fasta", arity: '1..*')

    exec:
    def sequences = uk.ac.ebi.interpro.FastaFile.chunkSequences(fasta, chunk_size, chunk_overlap)
    def fileIndex = 1
    def seqCount = 0
    def writer = task.workDir.resolve("tmbed_${fileIndex}.fasta").newWriter("UTF-8")
    sequences.each { seq ->
        if (batch_size > 0 && seqCount >= batch_size) {
            writer.close()
            fileIndex += 1
            seqCount = 0
            writer = task.workDir.resolve("tmbed_${fileIndex}.fasta").newWriter("UTF-8")
        }
        seq.chunks.eachWithIndex { chunk, idx ->
            writer.writeLine(">${seq.id}_${idx + 1}")
            chunk.eachMatch(/.{1,60}/) { writer.writeLine(it) }
            seqCount += 1
        }
    }

    writer.close()
}

process RUN_TMBED_CPU {
    label     'mem_high', 'time_medium', 'dynamic'
    container 'interpro/tmbed:1.0.2'

    input:
    tuple val(meta), path(fasta)
    val batch_size

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -t ${task.cpus} \
        -p tmbed.pred \
        --batch-size ${batch_size}
    """
}

process RUN_TMBED_GPU {
    label     'mem_medium', 'time_short', 'use_gpu'
    container 'interpro/tmbed:1.0.2'

    input:
    tuple val(meta), path(fasta)
    val batch_size

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -p tmbed.pred \
        --use-gpu \
        --batch-size ${batch_size}
    """
}

process PARSE_TMBED {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(tmbed_out)
    val chunk_overlap
    val smooth_window

    output:
    tuple val(meta), path("tmbed.json")

    exec:
    /**
     * Parse file into chunks grouped by base sequence ID
     * TMbed output format consists of repeating blocks of 3 lines, 
     * where each block corresponds to one sequence. The 3 lines are:
     * >seqId
     * sequence
     * prediction
     */
    def chunksbySeqId = [:].withDefault { [] }
    def seqId = null
    def chunkId = null
    def prediction = null
    def lineCounter = 0
    tmbed_out.eachLine { line ->
        line = line.trim()
        if (lineCounter % 3 == 0) {
            def m = line =~ /^>(.+)_([0-9]+)$/
            assert m.size() == 1
            seqId = m[0][1]
            chunkId = m[0][2] as int
        } else if (lineCounter % 3 == 2) {
            prediction = line
            chunksbySeqId[seqId] << [id: chunkId, pred: prediction]
            seqId = null
            chunkId = null
            prediction = null
        }

        lineCounter += 1
    }

    def mergedList = []
    def hits = [:].withDefault { [:] }
    def libRelease = new uk.ac.ebi.interpro.SignatureLibraryRelease("TMbed", "1.0.2")
    def MODEL_TYPES = [  // Sig(acc, name, desc, type, lib, entry)
        "b": new uk.ac.ebi.interpro.Signature("TMbeta_out-to-in", "Transmembrane beta strand (out-to-in)", null, "Region", libRelease, null),
        "B": new uk.ac.ebi.interpro.Signature("TMbeta_in-to-out" ,"Transmembrane beta strand (in-to-out)", null, "Region", libRelease, null),
        "h": new uk.ac.ebi.interpro.Signature("TMhelix_out-to-in", "Transmembrane alpha helix (out-to-in)", null, "Region", libRelease, null),
        "H": new uk.ac.ebi.interpro.Signature("TMhelix_in-to-out", "Transmembrane alpha helix (in-to-out)", null, "Region", libRelease, null),
        "S": new uk.ac.ebi.interpro.Signature("Signal_peptide", "Signal peptide", null, "Region", libRelease, null)
    ]
    chunksbySeqId.each { hitId, chunks ->
        // Sort chunks to ensure correct order before starting to merge 
        chunks = chunks.sort { a, b ->
            a.id <=> b.id
        }

        // Initialize merged prediction with the first chunk
        def mergedPred = new StringBuilder(chunks[0].pred)

        // Merge subsequent chunks
        if (chunks.size() > 1) {
            chunks.subList(1, chunks.size()).each { chunk ->
                def chunkPred = chunk.pred
                int chunkPredLen = chunkPred.length()

                // Effective overlap might be smaller than chunk_overlap for the last chunk
                int effOverlap = Math.min(chunk_overlap, chunkPredLen)
                int mergedLen = mergedPred.length()
                int overlapStart = mergedLen - effOverlap

                // Build a "reconciled" overlap region by comparing per-residue predictions
                def mergedOverlap = new StringBuilder(effOverlap)
                effOverlap.times { j ->
                    def a = mergedPred.charAt(overlapStart + j)
                    def b = chunkPred.charAt(j)
                    if (a == b) {
                        mergedOverlap.append(a)
                    } else {
                        /* Conflict between predictions
                           Resolve by local consensus using a +/- 2 residue window */
                        int winStart = Math.max(0, j - 2)
                        int winEnd = Math.min(effOverlap, j + 3)
                        def winPrev = mergedPred.substring(overlapStart + winStart, overlapStart + winEnd)
                        def winNext = chunkPred.substring(winStart, winEnd)
                        def mostPrev = mostCommonChar(winPrev)
                        def mostNext = mostCommonChar(winNext)
                        /* If the current prediction agrees with the local consensus
                           of the previous chunk, we prefer that prediction as
                           it reflects the accumulated consensus and preserves continuity
                           across chunk boundaries */
                        mergedOverlap.append(mostPrev == a ? mostPrev : mostNext)
                    }
                }

                // Replace the overlapping region in mergedPred with the reconciled overlap
                mergedPred.setLength(overlapStart)
                mergedPred.append(mergedOverlap).append(chunkPred, effOverlap, chunkPredLen)
            }
        }
        
        // Smoothing step to avoid short suprious regions introduced by merging
        def mergedPredStr = mergedPred.toString()
        if (chunks.size() > 1 && smooth_window > 1) {
            mergedPredStr = smoothString(mergedPredStr, smooth_window)
        }

        // Parse merged predictions into Match objects
        def currentMatch = null
        def start = null
        def end = null
        mergedPredStr.eachWithIndex { symbol, position ->
            end = position
            if (symbol != ".") { // Found a hit
                if (!currentMatch) {  // Start a new match
                    def signature = MODEL_TYPES[symbol]
                    currentMatch = hits[hitId].computeIfAbsent(signature.accession) {
                        def newMatch = new uk.ac.ebi.interpro.Match(signature.accession, signature)
                        newMatch
                    }
                    start = position + 1

                } else if (currentMatch.modelAccession != MODEL_TYPES[symbol].accession) { // Found a new hit
                    currentMatch.addLocation(new uk.ac.ebi.interpro.Location(start, position))
                    def signature = MODEL_TYPES[symbol]
                    currentMatch = hits[hitId].computeIfAbsent(signature.accession) {
                        def newMatch = new uk.ac.ebi.interpro.Match(signature.accession, signature)
                        newMatch
                    }
                    start = position + 1
                }
                // else, parsing another symbol from the same currentMatch/hit
            } else {
                if (currentMatch) {
                    currentMatch.addLocation(new uk.ac.ebi.interpro.Location(start, position))
                    currentMatch = null
                }
            }
        }

        // Add the final match for this sequence
        if (currentMatch) {
            currentMatch.addLocation(new uk.ac.ebi.interpro.Location(start, end))
        }
    }

    def filepath = task.workDir.resolve("tmbed.json")
    filepath.text  = groovy.json.JsonOutput.toJson(hits)
}


/** Apply categorical smoothing (mode over a sliding window) */
def smoothString(String s, int window) {
    int half = Math.floor(window / 2) as int
    return (0..<s.size()).collect { i ->
        int start = Math.max(0, i - half)
        int end = Math.min(s.size(), i + half + 1)
        mostCommonChar(s.substring(start, end))
    }.join()
}

/** Return the most frequent character in a string */
def mostCommonChar(String s) {
    def counts = [:].withDefault { 0 }
    s.each { counts[it] += 1 }
    return counts.max { it.value }.key
}
