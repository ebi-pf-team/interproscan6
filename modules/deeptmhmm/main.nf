process RUN_DEEPTMHMM_CPU {
    label       'mem_high'
    label       'time_medium'
    container   'interpro/deeptmhmm:1.0'
    stageInMode 'copy'

    input:
    tuple val(meta), path(fasta)
    path tmhmm_dir

    output:
    tuple val(meta), path("outdir")

    script:
    """
    # deeptmhmm has a hard coded assumption it is being run within its dir
    cd ${tmhmm_dir}
    python3 predict.py \
        --fasta ../${fasta} \
        --output-dir ../outdir
    cd ..
    rm -r ${tmhmm_dir} outdir/embeddings
    chmod -R 777 outdir
    """
}

process RUN_DEEPTMHMM_GPU {
    label       'mem_high'
    label       'time_short'
    label       'use_gpu'
    container   'interpro/deeptmhmm:1.0'
    stageInMode 'copy'

    input:
    tuple val(meta), path(fasta)
    path tmhmm_dir

    output:
    tuple val(meta), path("outdir")

    script:
    """
    # deeptmhmm has a hard coded assumption it is being run within its dir
    cd ${tmhmm_dir}
    python3 predict.py \
        --fasta ../${fasta} \
        --output-dir ../outdir
    cd ..
    rm -r ${tmhmm_dir} outdir/embeddings
    chmod -R 777 outdir
    """
}

process PARSE_DEEPTMHMM {
    label    'mem_low'
    label    'time_short'
    executor 'local'

    input:
    tuple val(meta), val(tmhmm_output)

    output:
    tuple val(meta), path("tmhmm.json")

    exec:
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("DeepTMHMM", "1.0")
    def MODEL_TYPES = [
        "Beta sheet": ["Transmembrane beta barrel", new uk.ac.ebi.interpro.Signature("Transmembrane beta barrel", null, null, "Region", library, null)],
        "periplasm": ["Periplasmic Domain", new uk.ac.ebi.interpro.Signature("Periplasmic Domain", null, null, "Region", library, null)],
        "signal": ["Signalp Peptide", new uk.ac.ebi.interpro.Signature("Signal Peptide", null, null, "Region", library, null)],
        "TMhelix": ["Transmembrane alpha helix", new uk.ac.ebi.interpro.Signature("Transmembrane alpha helix", null, null, "Region", library, null)],
    ]
    def hits = [:]
    def seqId
    file("${tmhmm_output}/TMRs.gff3").eachLine { line ->
        def lineData = line.split("\t")
        if (line.startsWith("//") || line.startsWith("#")) {  // stops '##gff-version 3' line breaking assert
            return
        }
        assert lineData.size() == 4
        if (MODEL_TYPES.containsKey(lineData[1])) { // e.g. tr_A0A009GMU8_A0A009GMU8_9GAMM periplasm 30 184
            seqId = lineData[0]
            hits.computeIfAbsent(seqId) { [:] }
            def (modelAcc, modelSig) = MODEL_TYPES[lineData[1]]
            hits[seqId].computeIfAbsent(modelAcc) {
                def match = new uk.ac.ebi.interpro.Match(modelAcc, modelSig)
                match
            }
            def start = lineData[2].toInteger()
            def end = lineData[3].toInteger()
            hits[seqId][modelAcc].addLocation(new uk.ac.ebi.interpro.Location(start, end))
        }
    }

    def filepath = task.workDir.resolve("tmhmm.json")
    filepath.text = groovy.json.JsonOutput.toJson(hits)
}
