process RUN_SIGNALP_CPU {
    label     'mem_medium'
    label     'time_medium'
    container 'interpro/signalp:6.0i'

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
        --torch_num_threads 1 \
        --write_procs 1 \
        --model_dir ${signalp_dir}/signalp-6-package/models
    rm -r signalp
    chmod -R 777 outdir
    """
}

process RUN_SIGNALP_GPU {
    label     'mem_medium'
    label     'time_medium'
    label     'use_gpu'
    container 'interpro/signalp:6.0i'

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
        --torch_num_threads 1 \
        --write_procs 1 \
        --model_dir ${signalp_dir}/signalp-6-package/models
    rm -r signalp
    chmod -R 777 outdir
    """
}

process PARSE_SIGNALP {
    label    'mem_low'
    label    'time_short'
    executor 'local'

    input:
    tuple val(meta), val(organism), val(mode), val(signalp_out)

    output:
    tuple val(meta), path("signalp.json")

    exec:
    def jsonPath = signalp_out.resolve("output.json")
    def jsonOutput = new groovy.json.JsonSlurper().parse(jsonPath)

    def modelAcc = "SignalP_${mode}_${organism}"
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("SignalP", "6.0i")
    def signatures = [
        "Sec/SPI"  : new uk.ac.ebi.interpro.Signature("SignalP-Sec-SPI", "Sec/SPI", "Sec signal peptide", null, "Region", library, null),
        "Sec/SPII" : new uk.ac.ebi.interpro.Signature("SignalP-Sec-SPII", "Sec/SPII", "Lipoprotein signal peptide", null, "Region", library, null),
        "Tat/SPI"  : new uk.ac.ebi.interpro.Signature("SignalP-Tat-SPI", "Tat/SPI", "Tat signal peptide", null, "Region", library, null),
        "Tat/SPII" : new uk.ac.ebi.interpro.Signature("SignalP-Tat-SPII", "Tat/SPII", "Tat lipoprotein signal peptide", null, "Region", library, null),
        "Sec/SPIII": new uk.ac.ebi.interpro.Signature("SignalP-Sec-SPIII", "Sec/SPIII", "Pilin signal peptide", null, "Region", library, null),
    ]

    def gff3Path = signalp_out.resolve("output.gff3")
    def hits = [:]
    gff3Path.eachLine { line ->
        if (line.startsWith("#")) {
            return
        }

        def fields = line.split(/\t/)
        assert fields.size() == 9
        def seqHeader = fields[0]
        def seqId = seqHeader.split(/\s+/)[0]
        def start = fields[3].toInteger()
        def end = fields[4].toInteger()
        if (start < 0 || end < 0) {
            return
        }
        def score = Double.parseDouble(fields[5])
        def prediction = jsonOutput["SEQUENCES"][seqHeader]["Prediction"]

        def matcher = prediction =~ /\(([a-zA-Z]+\/[a-zA-Z]+)\)/
        def spType = matcher.find() ? matcher.group(1) : null
        def signature = signatures.get(spType)
        assert signature != null
        def match = new uk.ac.ebi.interpro.Match(modelAcc, signature)
        def location = new uk.ac.ebi.interpro.Location(start, end)
        location.score = score
        match.addLocation(location)
        hits[seqId] = [(modelAcc) : match]
    }

    def outputFile = task.workDir.resolve("signalp.json")
    outputFile.text  = groovy.json.JsonOutput.toJson(hits)
}
