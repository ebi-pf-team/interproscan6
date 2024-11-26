import groovy.json.JsonOutput

process RUN_NCBIFAM {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.out")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -Z 61295632 --cut_tc \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmsearch.out
    """
}

process PARSE_NCBIFAM {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("ncbifam.json")    

    exec:
    def outputFilePath = task.workDir.resolve("ncbifam.json")
    def matches = HMMER3.parseOutput(hmmseach_out.toString())
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

