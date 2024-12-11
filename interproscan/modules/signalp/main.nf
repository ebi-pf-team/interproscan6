import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process RUN_SIGNALP {
    label 'medium', 'signalp_container'

    input:
    tuple val(meta), path(fasta)
    val organism
    val mode
    path signalp_dir

    output:
    tuple val(meta), val(organism), val(mode), path("outdir")

    script:
    """
    cp -Lr ${signalp_dir}/signalp-6-package/signalp signalp
    python -m signalp.predict \
        --fastafile ${fasta} \
        --output_dir outdir \
        --format none \
        --organism ${organism} \
        --mode ${mode} \
        --write_procs 1 \
        --model_dir ${signalp_dir}/signalp-6-package/models
    rm -r signalp
    chmod -R 777 outdir
    """
}

process PARSE_SIGNALP {
    label 'small'

    input:
    tuple val(meta), val(organism), val(mode), val(signalp_out)

    output:
    tuple val(meta), path("signalp.json")

    exec:
    String signalDir = signalp_out.toString()

    File jsonFile = new File(signalDir, "output.json")
    def jsonSlurper = new JsonSlurper()
    def jsonOutput = jsonSlurper.parse(jsonFile)

    String modelAcc = "SignalP_${mode}_${organism}"

    def hits = [:]
    new File(signalDir, "output.gff3").eachLine { line ->
        if (line.startsWith("#")) {
            return
        }

        def fields = line.split(/\t/)
        assert fields.size() == 9
        String seqHeader = fields[0]
        String seqId = seqHeader.split(/\s+/)[0]
        int start = fields[3].toInteger()
        int end = fields[4].toInteger()
        Double score = Double.parseDouble(fields[5])
        String prediction = jsonOutput["SEQUENCES"][seqHeader]["Prediction"]

        SignatureLibraryRelease library = new SignatureLibraryRelease("SignalP", "6.0h")
        Match match = new Match(modelAcc)       
        match.signature = new Signature("SignalP", library)
        Location location = new Location(start, end, prediction)
        location.score = score
        match.addLocation(location)
        hits[seqId] = [(modelAcc) : match]
    }

    def outputFilePath = task.workDir.resolve("signalp.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
