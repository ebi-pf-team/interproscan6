import groovy.json.JsonSlurper
import groovy.json.JsonOutput

workflow AGGREGATE_RESULTS {
    take:
    all_results

    main:
    ch_aggregated_results = all_results
        .collect()
        .map { jsonFiles ->
            def combined = [:]
            jsonFiles.each { jsonFile ->
                def file = new File(jsonFile.toString())
                if (!file.exists()) {
                    println "File not found: ${jsonFile}"
                    log.error ("File not found: ${jsonFile}")
                    exit 5
                }
                try {
                    def json = new groovy.json.JsonSlurper().parse(file)
                    combined += json
                } catch (groovy.json.JsonException e) {
                    println "Corrupted or improperly formatted JSON file: ${jsonFile}"
                    log.error ("Corrupted or improperly formatted JSON file: ${jsonFile}", e)
                    exit 22
                }
            }
            def aggregated_result = JsonOutput.toJson(combined)
            def tempFile = File.createTempFile("aggregated_result", ".json")
            tempFile.write(JsonOutput.prettyPrint(aggregated_result))
            return tempFile.path
        }

    emit:
    ch_aggregated_results
}
