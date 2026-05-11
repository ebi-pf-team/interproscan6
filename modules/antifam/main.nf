import groovy.json.JsonOutput
import uk.ac.ebi.interpro.HMMER3

process PARSE_ANTIFAM {
    label    'mem_min','time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("antifam.json")

    exec:
    def filepath = task.workDir.resolve("antifam.json")
    def matches = HMMER3.parseOutput(hmmseach_out, "AntiFam")
    filepath.text = JsonOutput.toJson(matches)
}
