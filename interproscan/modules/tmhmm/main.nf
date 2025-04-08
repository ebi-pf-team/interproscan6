import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process RUN_DEEPTMHMM {
    label 'medium', 'deeptmhmm_container'
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
    label 'run_locally', 'deeptmhmm_container'

    input:
    tuple val(meta), val(tmhmm_output)

    output:
    tuple val(meta), path("tmhmm.json")

    exec:
    SignatureLibraryRelease library = new SignatureLibraryRelease("DeepTMHMM", "1.0")
    def MODEL_TYPES = [
        "Beta sheet": ["Transmembrane beta barrel", new Signature("Transmembrane beta barrel", library)],
        "periplasm": ["Periplasmic Domain", new Signature("Periplasmic Domain", library)],
        "signal": ["Signalp Peptide", new Signature("Signal Peptide", library)],
        "TMhelix": ["Transmembrane alpha helix", new Signature("Transmembrane alpha helix", library)],
    ]
    String tmhmmDir = tmhmm_output.toString()
    Map<String, Match> hits = [:]
    String seqId
    file("${tmhmm_output}/TMRs.gff3").eachLine { line ->
        def lineData = line.split("\t")
        if (line.startsWith("//") || line.startsWith("#")) {  // stops '##gff-version 3' line breaking assert
            return
        }
        assert lineData.size() == 4
        if (MODEL_TYPES.containsKey(lineData[1])) { // e.g. tr_A0A009GMU8_A0A009GMU8_9GAMM periplasm 30 184
            seqId = lineData[0]
            hits.computeIfAbsent(seqId) { [:] }
            (modelAcc, modelSig) = MODEL_TYPES[lineData[1]]
            hits[seqId].computeIfAbsent(modelAcc) {
                Match match = new Match(modelAcc)
                match.signature = modelSig
                match
            }
            int start = lineData[2].toInteger()
            int end = lineData[3].toInteger()
            hits[seqId][modelAcc].addLocation(new Location(start, end))
        }
    }

    def outputFilePath = task.workDir.resolve("tmhmm.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
