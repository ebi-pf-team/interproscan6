import groovy.json.JsonOutput

process RUN_NCBIFAM {
    label 'small'

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
    label 'small'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("ncbifam.json")

    exec:
    def outputFilePath = task.workDir.resolve("ncbifam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString())

    def processedMatches = hmmerMatches.collectEntries { seqId, matches ->
        [seqId, matches.collectEntries { modelAccession, match ->
            def updatedModelAccession = modelAccession.split("\\.")[0]
            match.modelAccession = updatedModelAccession
            [(updatedModelAccession): match]
        }]
    }

    def json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}

