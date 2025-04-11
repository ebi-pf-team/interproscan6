import groovy.json.JsonOutput

process PARSE_ANTIFAM {
    label 'run_locally'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("antifam.json")

    exec:
    def outputFilePath = task.workDir.resolve("antifam.json")
    def matches = HMMER3.parseOutput(hmmseach_out.toString(), "AntiFam")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
