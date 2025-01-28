import groovy.json.JsonOutput

process RUN_ANTIFAM {
    label 'small', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        --cut_ga \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process PARSE_ANTIFAM {
    label 'local'

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

