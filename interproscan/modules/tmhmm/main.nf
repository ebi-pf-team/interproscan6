import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process RUN_DEEPTMHMM {
    label 'medium', 'deeptmhmm_runner'
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
    rm -r ${tmhmm_dir}
    chmod -R 777 outdir
    """
}

process PARSE_DEEPTMHMM {
    label 'deeptmhmm_runner'

    input:
    tuple val(meta), val(tmhmm_output)

    output:
    tuple val(meta), path("tmhmm.json")

    exec:
    Map<String, String> MODEL_TYPES = [
        "signal": "Signalp Peptide",
        "TMhelix": "Transmembrane Helix",
        "periplasm": "Periplasmic Domain",
        "Beta": "Beta Sheet"
    ]
    SignatureLibraryRelease library = new SignatureLibraryRelease("DeepTMHMM", "1.0")
    Signature betaSheetSig = new Signature("Beta Sheet", library)
    Signature periplasmSig = new Signature("Periplasmic Domain", library)
    Signature signalPeptideSig = new Signature("Signal Peptide", library)
    Signature tmHelixSig = new Signature("Transmembrane Helix", library)

    String tmhmmDir = tmhmm_output.toString()
    Map<String, Match> hits = [:]
    String seqId
    file("${tmhmm_output}/TMRs.gff3").eachLine { line ->
        def lineData = line.split("\\s+")
        if (line.startsWith("//") || line.startsWith("##")) {  // MODEL_TYPES.con(lineData[1]) will crash on these linse
            return
        } else if (line.startsWith("# ")) { // e.g. # tr_A0A009GMU8_A0A009GMU8_9GAMM Number of predicted TMRs: 22
            seqId = lineData[1]
        } else if (MODEL_TYPES.containsKey(lineData[1])) { // e.g. tr_A0A009GMU8_A0A009GMU8_9GAMM periplasm 30 184
            hits.computeIfAbsent(seqId) { [:] }
            String modelAcc = MODEL_TYPES[lineData[1]]
            hits[seqId].computeIfAbsent(modelAcc) {
                new Match(modelAcc).with {
                    signature = switch (modelAcc) {
                        case "Beta Sheet" -> betaSheetSig
                        case "Periplasmic Domain" -> periplasmSig
                        case "Signal Peptide" -> signalPeptideSig
                        case "Transmembrane Helix" -> tmHelixSig
                    }
                    it
                }
            }
            int start = lineData[-2].toInteger()
            int end = lineData[-1].toInteger()
            hits[seqId][modelAcc].addLocation(new Location(start, end))
        }
    }

    def outputFilePath = task.workDir.resolve("tmhmm.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
