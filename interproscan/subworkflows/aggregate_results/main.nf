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
                def json = new groovy.json.JsonSlurper().parse(new File(jsonFile.toString()))
                combined += json
            }
            def aggregated_result = JsonOutput.toJson(combined)
            def tempFile = File.createTempFile("aggregated_result", ".json")
            tempFile.write(JsonOutput.prettyPrint(aggregated_result))
            return tempFile.path
        }

    emit:
    ch_aggregated_results
}
