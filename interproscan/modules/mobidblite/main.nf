import groovy.json.JsonOutput

process RUN_MOBIDBLITE {
    label 'mobidblite_runner'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("output.tsv")

    script:
    """
    idrpred ${fasta} output.tsv
    """
}


process PARSE_MOBIDBLITE {
    label 'analysis_parser'

    input:
    tuple val(meta), val(mobidblite_output)

    output:
    tuple val(meta), path("mobidblite.json")

    exec:
    def outputFilePath = task.workDir.resolve("mobidblite.json")
    def matches = [:]

    file(mobidblite_output.toString()).eachLine { line ->
        def lineData = line.split(/\s+/)
        def sequenceId = lineData[0]
        def start = lineData[1].toInteger()
        def end = lineData[2].toInteger()
        def feature = lineData[3] ?: ""

        if (matches.containsKey(sequenceId)) {
            matches[sequenceId]["mobidb_lite"]["locations"].add([
                start: start,
                end: end,
                "sequence-feature": feature
            ])
        } else {
            matches[sequenceId] = [
                "mobidb_lite": [
                    member_db: "mobidb_lite",
                    accession: "mobidb-lite",
                    name: "disorder_prediction",
                    description: "consensus disorder prediction",
                    locations: [[
                        start: start,
                        end: end,
                        "sequence-feature": feature
                    ]]
                ]
            ]
        }
    }

    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
