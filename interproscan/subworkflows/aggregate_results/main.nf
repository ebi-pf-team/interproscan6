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
                try {
                    def json = new groovy.json.JsonSlurper().parse(new File(jsonFile.toString()))
                    combined += json
                } catch (FileNotFoundException e) {
                    println "File not found: ${jsonFile}"
                    log.error ("File not found: ${jsonFile}", e)
                    exit 5
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
