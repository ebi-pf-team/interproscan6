import groovy.json.JsonOutput

process RUN_SIGNALP {
    label 'signalp_runner'

    input:
    tuple val(meta), path(fasta)
    val orgType
    val mode

    output:
    tuple val(meta), path("signalp_out"), val(orgType)

    script:
    """
    signalp6 \
        --organism ${orgType} \
        --fastafile ${fasta} \
        --output_dir signalp_out \
        --mode ${mode}
    """
}

process PARSE_SIGNALP {
    label 'analysis_parser'

    input:
    tuple val(meta), val(signalp_out), val(orgType)
    val threshold

    output:
    path "signalp.json"

    exec:
    def outputFilePath = task.workDir.resolve("signalp.json")
    def matches = SIGNALP.parseOutput(signalp_out.toString(), threshold, orgType)
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
