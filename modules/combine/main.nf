import com.fasterxml.jackson.databind.ObjectMapper
import groovy.json.JsonOutput

process COMBINE_MATCHES_LOCAL {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path("combined_matches.json")

    exec:
    def aggregatedMatches = [:]  // seqMd5: {modelAcc: Match} -- stops a seqMd5 appearing multiple times in the output
    members_matches.each { matchesPath ->
        if (matchesPath) {  // can be null if LOOKUP_MATCHES failed or was skipped
            new ObjectMapper().readValue(new File(matchesPath.toString()), Map).each { seqMd5, matches ->
                aggregatedMatches.computeIfAbsent(seqMd5, { [:] }).putAll(matches)
            }
        }
    }

    String outputFilePath = task.workDir.resolve("combined_matches.json")
    def json = JsonOutput.toJson(aggregatedMatches)
    new File(outputFilePath.toString()).write(json)
}

process COMBINE_MATCHES {
    label    'mem_high', 'time_medium'

    input:
    tuple val(meta), val(members_matches)

    output:
    tuple val(meta), path("combined_matches.json")

    script:
    """
    #!/usr/bin/env python

    import glob
    import json
    
    files = glob.glob("*.json")
    results = {}
    for file in files:
         if file:
            with open(file, "rt") as fh:
                matches = json.load(fh)

            for seq_id in matches:
                try:
                    obj = results[seq_id]
                except KeyError:
                    obj = results[seq_id] = {}

                obj.update(matches[seq_id])

    with open("combined_matches.json", "wt") as fh:
        json.dump(results, fh)
    """
}