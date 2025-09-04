import groovy.json.JsonOutput

import Match

process RUN_TMBED_CPU {
    label 'large', 'dynamic', 'tmbed_container'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -t ${task.cpus} \
        -p tmbed.pred
    """
}

process RUN_TMBED_GPU {
    label 'large', 'tmbed_container', 'use_gpu'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("tmbed.pred")

    script:
    """
    tmbed predict \
        -f ${fasta} \
        -m /opt/tmbed/models \
        -p tmbed.pred \
        --use-gpu
    """
}

process PARSE_TMBED {
    label 'tiny'
    executor 'local'

    input:
    tuple val(meta), val(tmbed_out)

    output:
    tuple val(meta), path("tmbed.json")

    exec:
    SignatureLibraryRelease LIBRARY = new SignatureLibraryRelease("TMbed", "1.0.2")
    def MODEL_TYPES = [  // Sig(acc, name, desc, type, lib, entry)
        "b": new Signature("TMbeta_out-to-in", "Transmembrane beta strand (out-to-in)", null, "Region", LIBRARY, null),
        "B": new Signature("TMbeta_in-to-out" ,"Transmembrane beta strand (in-to-out)", null, "Region", LIBRARY, null),
        "h": new Signature("TMhelix_out-to-in", "Transmembrane alpha helix (out-to-in)", null, "Region", LIBRARY, null),
        "H": new Signature("TMhelix_in-to-out", "Transmembrane alpha helix (in-to-out)", null, "Region", LIBRARY, null),
        "S": new Signature("Signal_peptide", "Signal peptide", null, "Region", LIBRARY, null)
    ]


    Map<String, Map<String, Match>> hits = [:]
    Match currentMatch = null
    def seqMd5 = null
    def start = null
    def end = null
    def lineCounter = 0
    /* TMbed output:
    >seqId
    sequence
    prediction */
    new File(tmbed_out.toString()).eachLine { line ->
        line = line.trim()
        if (lineCounter % 3 == 0) {
            seqMd5       = line.substring(1)
            currentMatch = null
            hits.computeIfAbsent(seqMd5) { [:] }
        } else if (lineCounter % 3 == 2) {
            line.eachWithIndex { symbol, position ->
                end = position
                if (symbol != ".") { // Found a hit
                    if (!currentMatch) {  // Start a new match
                        Signature signature = MODEL_TYPES[symbol]
                        currentMatch = hits[seqMd5].computeIfAbsent(signature.accession) {
                            Match newMatch = new Match(signature.accession)
                            newMatch.signature = signature
                            newMatch
                        }
                        start = position + 1

                    } else if (currentMatch.modelAccession != MODEL_TYPES[symbol].accession) { // Found a new hit
                        currentMatch.addLocation(new Location(start, position))
                        Signature signature = MODEL_TYPES[symbol]
                        currentMatch = hits[seqMd5].computeIfAbsent(signature.accession) {
                            Match newMatch = new Match(signature.accession)
                            newMatch.signature = signature
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
        lineCounter += 1
    }

    def outputFilePath = task.workDir.resolve("tmbed.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
