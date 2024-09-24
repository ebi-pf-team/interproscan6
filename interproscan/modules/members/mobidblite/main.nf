import groovy.json.JsonOutput

process MOBIDBLITE_RUNNER {
    label 'mobidblite_runner'

    input:
    path "fasta"

    output:
    path "mobidb_out.tsv"

    script:
    """
    idrpred ${fasta} mobidb_out.tsv
    """
}


process MOBIDBLITE_PARSER {
    label 'analysis_parser'

    input:
    val mobidb_out

    output:
    path "mobidb_parsed.json"

    exec:
    def outputFilePath = task.workDir.resolve("mobidb_parsed.json")
    def matches = [:]

    file(mobidb_out.toString()).eachLine { line ->
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

    new File(outputFilePath.toString()).text = groovy.json.JsonOutput.toJson(matches)
}
