import groovy.json.JsonOutput


process RUN_PRINTS {
    label 'prints_runner'

    input:
    tuple val(meta), path(fasta)
    val pvalFile

    output:
    tuple val(meta), path("prints_output")

    script:
    """
    $projectDir/bin/prints/fingerPRINTScan \
        ${pvalFile} \
        ${fasta} \
        -e 0.0001 -d 10 -E 257043 84355444 -fj -o 15 > prints_output
    """
}

process PARSE_PRINTS {
    label 'analysis_parser'

    input:
    tuple val(meta), path(prints_output)
    val hierarchy

    output:
    tuple val(meta), path("prints.json")

    exec:
    def outputFilePath = task.workDir.resolve("prints.json")
    def matches = PRINTS.parseOutput(prints_output.toString(), hierarchy)
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}