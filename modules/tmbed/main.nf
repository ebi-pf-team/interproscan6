// codenarc-disable AllowedDirectivesRule
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease

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
    def sequences = FastaFile.chunkSequences(fasta, chunk_size, chunk_overlap)
    def fileIndex = 1
    def seqCount = 0
    def writer = null
    def openNewFile = {
        if (writer) writer.close()
        writer = task.workDir.resolve("tmbed_${fileIndex}.fasta").newWriter("UTF-8")
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
            writer.writeLine(">${seq.id}_${idx + 1}")
            chunk.eachMatch(/.{1,60}/) { writer.writeLine(it) }
            seqCount++
        }
    }

    if (writer) writer.close()
}

process RUN_TMBED_CPU {
    label 'mem_high', 'time_medium', 'dynamic', 'tmbed_container'

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
    label 'mem_medium', 'time_short', 'tmbed_container', 'use_gpu'

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
    SignatureLibraryRelease libRelease = new SignatureLibraryRelease("TMbed", "1.0.2")
    def MODEL_TYPES = [  // Sig(acc, name, desc, type, lib, entry)
        "b": new Signature("TMbeta_out-to-in", "Transmembrane beta strand (out-to-in)", null, "Region", libRelease, null),
        "B": new Signature("TMbeta_in-to-out" ,"Transmembrane beta strand (in-to-out)", null, "Region", libRelease, null),
        "h": new Signature("TMhelix_out-to-in", "Transmembrane alpha helix (out-to-in)", null, "Region", libRelease, null),
        "H": new Signature("TMhelix_in-to-out", "Transmembrane alpha helix (in-to-out)", null, "Region", libRelease, null),
        "S": new Signature("Signal_peptide", "Signal peptide", null, "Region", libRelease, null)
    ]
    chunksbySeqId.each { hitId, chunks ->
        // Sort chunks to ensure correct order before starting to merge 
        chunks = chunks.sort { a, b ->
            a.id <=> b.id
        }

        // Initialize merged prediction with the first chunk
        def mergedPred = new StringBuilder(chunks[0].pred)

        // Merge subsequent chunks
        for (int i = 1; i < chunks.size(); i++) {
            def chunkPred = chunks[i].pred

            // Effective overlap might be smaller than chunk_overlap for the last chunk
            int effOverlap = Math.min(chunk_overlap, chunkPred.length())

            // Extract overlapping regions
            def overlapPrev = mergedPred.substring(mergedPred.length() - effOverlap)
            def overlapNext = chunkPred.substring(0, effOverlap)

            // Build a "reconciled" overlap region by comparing per-residue predictions
            def mergedOverlap = new StringBuilder()
            for (int j = 0; j < effOverlap; j++) {
                def a = overlapPrev.charAt(j)
                def b = overlapNext.charAt(j)
                if (a == b) {
                    mergedOverlap.append(a)
                } else {
                    /* Conflict between predictions
                       Resolve by local consensus using a +/- 2 residue window */
                    def winPrev = overlapPrev.substring(Math.max(0, j - 2), Math.min(effOverlap, j + 3))
                    def winNext = overlapNext.substring(Math.max(0, j - 2), Math.min(effOverlap, j + 3))
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
            mergedPred.setLength(mergedPred.length() - effOverlap)
            mergedPred.append(mergedOverlap).append(chunkPred.substring(effOverlap))
        }

        // Smoothing step to avoid short suprious regions introduced by merging
        def mergedPredStr = mergedPred.toString()
        if (chunks.size() > 1 && smooth_window > 1) {
            mergedPredStr = smoothString(mergedPredStr, smooth_window)
        }

        // Parse merged predictions into Match objects
        Match currentMatch = null
        def start = null
        def end = null
        mergedPredStr.eachWithIndex { symbol, position ->
            end = position
            if (symbol != ".") { // Found a hit
                if (!currentMatch) {  // Start a new match
                    Signature signature = MODEL_TYPES[symbol]
                    currentMatch = hits[hitId].computeIfAbsent(signature.accession) {
                        Match newMatch = new Match(signature.accession, signature)
                        newMatch
                    }
                    start = position + 1

                } else if (currentMatch.modelAccession != MODEL_TYPES[symbol].accession) { // Found a new hit
                    currentMatch.addLocation(new Location(start, position))
                    Signature signature = MODEL_TYPES[symbol]
                    currentMatch = hits[hitId].computeIfAbsent(signature.accession) {
                        Match newMatch = new Match(signature.accession, signature)
                        newMatch
                    }
                    start = position + 1
                }
                // else, parsing another symbol from the same currentMatch/hit
            } else {
                if (currentMatch) {
                    currentMatch.addLocation(new Location(start, position))
                    currentMatch = null
                }
            }
        }

        // Add the final match for this sequence
        if (currentMatch) {
            currentMatch.addLocation(new Location(start, end))
        }
    }

    def filepath = task.workDir.resolve("tmbed.json")
    filepath.text = JsonOutput.toJson(hits)
}


/** Apply categorical smoothing (mode over a sliding window) */
def smoothString(String s, int window) {
    int half = Math.floor(window / 2) as int
    def out = new StringBuilder()
    for (int i = 0; i < s.size(); i++) {
        int start = Math.max(0, i - half)
        int end = Math.min(s.size(), i + half + 1)
        def windowSeq = s.substring(start, end)
        out.append(mostCommonChar(windowSeq))
    }
    return out.toString()
}

/** Return the most frequent character in a string */
def mostCommonChar(String s) {
    def counts = [:].withDefault { 0 }
    s.each { counts[it]++ }
    return counts.max { it.value }.key
}