import groovy.json.JsonOutput

process PARSE_NCBIFAM {
    label 'native'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("ncbifam.json")

    exec:
    def outputFilePath = task.workDir.resolve("ncbifam.json")
    def hmmerMatches = HMMER3.parseOutput(hmmseach_out.toString(), "NCBIFAM")

    def processedMatches = hmmerMatches.collectEntries { seqId, matches ->
        [seqId, matches.collectEntries { modelAccession, match ->
            def updatedModelAccession = modelAccession.split("\\.")[0]
            match.modelAccession = updatedModelAccession
            match.signature.accession = updatedModelAccession
            [(updatedModelAccession): match]
        }]
    }

    def json = JsonOutput.toJson(processedMatches)
    new File(outputFilePath.toString()).write(json)
}
