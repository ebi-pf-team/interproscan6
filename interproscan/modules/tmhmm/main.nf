import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process RUN_TMHMM {
    label 'medium', 'tmhmm_runner'
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
    chmod -R 777 outdir
    """
}

process PARSE_TMHMM {
    label 'analysis_parser'

    input:
    tuple val(meta), val(tmhmm_output)

    output:
    tuple val(meta), path("tmhmm.json")

    exec:
    List<String> NON_TMRS = ["inside", "outside"]
    Map<String, String> MODEL_TYPES = [
        "signal": "SIGNAL_PEPTIDE",
        "TMhelix": "TRANSMEMBRANE_HELIX",
        "periplasm": "PERIPLASMIC_DOMAIN",
        "Beta": "BETA_SHEET"
    ]

    String tmhmmDir = tmhmm_output.toString()
    Map<String, Match> hits = [:]
    String seqId
    file("${tmhmm_output}/TMRs.gff3").eachLine { line ->
        def lineData = line.split("\\s+")
        if (line.startsWith("//") || line.startsWith("##") || (line.startsWith("#") && line ==~ /#.*Length: \d+$/)) {
            return
        } else if (line ==~ ~/^#.*TMRs: \d+$/) { // e.g. # tr_A0A009GMU8_A0A009GMU8_9GAMM Number of predicted TMRs: 22
           seqId = lineData[1]
        } else if (!NON_TMRS.contains(lineData[1])) { // e.g. tr_A0A009GMU8_A0A009GMU8_9GAMM periplasm 30 184
            hits.computeIfAbsent(seqId) { [:] }
            String modelAcc = MODEL_TYPES.get(lineData[1], lineData[1].toUpperCase().replace(" ", "_"))
            hits[seqId].computeIfAbsent(modelAcc) {
                Match match = new Match(modelAcc)
                SignatureLibraryRelease library = new SignatureLibraryRelease("tmhmm", "1.0")
                match.signature = new Signature(modelAcc, library)
                match
            }
            int start = lineData[-2].toInteger()
            int end = lineData[-1].toInteger()
            Location location = new Location(start, end)
            location.fragments.add(new LocationFragment(start, end, "CONTINUOUS"))
            hits[seqId][modelAcc].addLocation(location)
        }
    }

    def outputFilePath = task.workDir.resolve("tmhmm.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
