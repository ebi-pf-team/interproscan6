import groovy.json.JsonSlurper
import groovy.json.JsonOutput

workflow AGGREGATE_RESULTS {
    take:
    all_results

    main:
    ch_aggregated_results = all_results
        .collect()
        .map { jsonFiles ->
            def merged = [:]
            jsonFiles.each { jsonFile ->
                def file = new File(jsonFile.toString())
                if (!file.exists()) {
                    log.error ("File not found: ${jsonFile}")
                    exit 5
                }
                try {
                    def hits = new groovy.json.JsonSlurper().parse(file)
                    hits.each { key, value ->
                        merged[key] = merged.containsKey(key) ? merged[key] + value : value
                    }
                } catch (groovy.json.JsonException e) {
                    log.error ("Corrupted or improperly formatted JSON file: ${jsonFile}", e)
                    exit 22
                }
            }
            def aggregated_result = JsonOutput.toJson(merged)
            def outputFile = new File("${workDir}/aggregated_result.json")
            outputFile.write(aggregated_result)
            return outputFile.path
        }

    emit:
    ch_aggregated_results
}
