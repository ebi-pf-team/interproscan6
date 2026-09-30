process RUN_DEEPTMHMM2_CPU {
    label     'mem_high'
    label     'time_medium'
    label     'dynamic'
    container 'interpro/deeptmhmm2:0.1.0-ffccb76'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("deeptmhmm2.json")

    script:
    """
    dtm2 \
        --model-dir /opt/deeptmhmm2 \
        --device cpu \
        --threads ${task.cpus} \
        --simplify-io \
        ${fasta} \
        outdir

    mv outdir/predictions.json deeptmhmm2.json
    rm -r outdir
    """
}

process RUN_DEEPTMHMM2_GPU {
    label     'mem_high'
    label     'time_short'
    label     'use_gpu'
    container 'interpro/deeptmhmm2:0.1.0-ffccb76'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("deeptmhmm2.json")

    script:
    """
    dtm2 \
        --model-dir /opt/deeptmhmm2 \
        --device cuda \
        --simplify-io \
        ${fasta} \
        outdir

    mv outdir/predictions.json deeptmhmm2.json
    rm -r outdir
    """
}

process PARSE_DEEPTMHMM2 {
    label    'mem_low'
    label    'time_short'
    executor 'local'

    input:
    tuple val(meta), val(dtm2_out)

    output:
    tuple val(meta), path("deeptmhmm2.json")

    exec:
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("DeepTMHMM2", "0.1.0")
    /* Segment names as reported with --simplify-io.
       Inside/outside follow the secretory-pathway convention: "outside" is the side
       the signal peptide is translocated to (e.g. extracellular or lumenal) */
    def SEGMENTS = [  // Sig(acc, name, desc, type, lib, entry)
        "signal"         : new uk.ac.ebi.interpro.Signature("SIGNAL_PEPTIDE", "Signal peptide", null, "Region", library, null),
        "transit_peptide": new uk.ac.ebi.interpro.Signature("TRANSIT_PEPTIDE", "Transit peptide", null, "Region", library, null),
        "TMhelix"        : new uk.ac.ebi.interpro.Signature("TM_HELIX", "Transmembrane alpha helix", null, "Region", library, null),
        "Beta sheet"     : new uk.ac.ebi.interpro.Signature("TM_BETA", "Transmembrane beta strand", null, "Region", library, null),
        "reentrant"      : new uk.ac.ebi.interpro.Signature("REENTRANT", "Re-entrant region", null, "Region", library, null),
        "interfacial"    : new uk.ac.ebi.interpro.Signature("INTERFACIAL", "Interfacial helix", null, "Region", library, null),
        "inside"         : new uk.ac.ebi.interpro.Signature("INSIDE", "Inside",
                                    "Non-membrane region on the side opposite to the signal peptide (e.g. cytoplasmic for plasma membrane proteins)", "Region", library, null),
        "outside"        : new uk.ac.ebi.interpro.Signature("OUTSIDE", "Outside",
                                    "Non-membrane region on the side of the signal peptide (e.g. extracellular or lumenal)", "Region", library, null),
    ]
    def TOPOLOGY_SIDES = ["inside", "outside"]

    def hits = [:]
    def predictions = new groovy.json.JsonSlurper().parseText(dtm2_out.text)
    predictions.each { prediction ->
        if (prediction.containsKey("metadata")) {
            return
        }

        // Only report inside/outside regions for membrane proteins
        def is_membrane = !prediction.type.startsWith("Globular")
        def matches = [:]
        prediction.segments.each { segment ->
            def (name, start, end) = segment
            if (!SEGMENTS.containsKey(name)) {
                throw new Exception("Unknown DeepTMHMM2 segment '${name}' for sequence ${prediction.id}")
            } else if (!is_membrane && TOPOLOGY_SIDES.contains(name)) {
                return
            }

            def signature = SEGMENTS[name]
            def match = matches.computeIfAbsent(signature.accession) {
                new uk.ac.ebi.interpro.Match(signature.accession, signature)
            }
            match.addLocation(new uk.ac.ebi.interpro.Location(start as int, end as int))
        }

        if (matches) {
            hits[prediction.id] = matches
        }
    }

    def filepath = task.workDir.resolve("deeptmhmm2.json")
    filepath.text = groovy.json.JsonOutput.toJson(hits)
}
