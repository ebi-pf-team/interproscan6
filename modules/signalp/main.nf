import groovy.json.JsonSlurper
import groovy.json.JsonOutput

import Match

process RUN_SIGNALP_CPU {
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
    sed -i "s|'-tt',\\s*|'-tt', type=int, |" signalp/predict.py
    python -m signalp.predict \
        --fastafile ${fasta} \
        --output_dir outdir \
        --format none \
        --organism ${organism} \
        --mode ${mode} \
        --torch_num_threads 1 \
        --write_procs 1 \
        --model_dir ${signalp_dir}/signalp-6-package/models
    rm -r signalp
    chmod -R 777 outdir
    """
}

process RUN_SIGNALP_GPU {
    label 'medium', 'signalp_container', 'use_gpu'

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
    sed -i "s|'-tt',\\s*|'-tt', type=int, |" signalp/predict.py
    python -m signalp.predict \
        --fastafile ${fasta} \
        --output_dir outdir \
        --format none \
        --organism ${organism} \
        --mode ${mode} \
        --torch_num_threads 1 \
        --write_procs 1 \
        --model_dir ${signalp_dir}/signalp-6-package/models
    rm -r signalp
    chmod -R 777 outdir
    """
}

process PARSE_SIGNALP {
    label    'tiny'
    executor 'local'

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
    SignatureLibraryRelease library = new SignatureLibraryRelease("SignalP", "6.0h")
    def signatures = [
        "Sec/SPI"  : new Signature("SignalP-Sec-SPI", "Sec/SPI", "Sec signal peptide", library, null),
        "Sec/SPII" : new Signature("SignalP-Sec-SPII", "Sec/SPII", "Lipoprotein signal peptide", library, null),
        "Tat/SPI"  : new Signature("SignalP-Tat-SPI", "Tat/SPI", "Tat signal peptide", library, null),
        "Tat/SPII" : new Signature("SignalP-Tat-SPII", "Tat/SPII", "Tat lipoprotein signal peptide", library, null),
        "Sec/SPIII": new Signature("SignalP-Sec-SPIII", "Sec/SPIII", "Pilin signal peptide", library, null),
    ]

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

        def matcher = prediction =~ /\(([a-zA-Z]+\/[a-zA-Z]+)\)/
        String spType = matcher.find() ? matcher.group(1) : null
        Signature signature = signatures.get(spType)
        assert signature != null
        Match match = new Match(modelAcc)    
        match.signature = signature
        Location location = new Location(start, end)
        location.score = score
        match.addLocation(location)
        hits[seqId] = [(modelAcc) : match]
    }

    def outputFilePath = task.workDir.resolve("signalp.json")
    def json = JsonOutput.toJson(hits)
    new File(outputFilePath.toString()).write(json)
}
