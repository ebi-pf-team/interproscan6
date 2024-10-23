import groovy.io.FileType

process PFSEARCH_RUNNER {
    label 'prosite_pfsearch_runner'

    input:
        tuple val(meta), path(fasta)
        val models_dir
        val skip_flagged_profiles

    output:
        path "prosite_profiles.out"
        path skip_flagged_profiles

    exec:
    // Run pfsearchV3 for all PROSITE profiles
    def outputFilePath = task.workDir.resolve("pfsearch.json")
    def outputFile = new File(outputFilePath.toString())

    List<String> profilePaths = []
    new File(models_dir).eachFileRecurse(FileType.FILES) { file ->
        profilePaths << file.absolutePath
    }
    profilePaths.each { profile ->
        def runCmd = ["/opt/pftools/pfsearchV3", profile, fasta, "-f", "-o", "7", "-t", "4"]
        try {
            def process = runCmd.execute()
            def output = process.text
            if (output.trim()) {
                outputFile.append(output + '\n')
            }
        } catch (Exception e) {
            println "Error running pfsearchV3 using cmd: ${runCmd.join(' ')}\nError: ${e.message}"
            System.exit(1)
        }
    }
}


process PFSEARCH_PARSER {
    label 'analysis_parser'

    input:
        tuple val(meta), path(pfsearch_out)
        path skip_flagged_profiles

    output:
        tuple val(meta), path("${pfsearch_out}-filtered.json")

    exec:
    Map profilesMatches = [:]

    with open(skip_flagged_profiles, 'r') as fh:
        to_skip = fh.read().splitlines()
    return to_skip

    def outputFilePath = task.workDir.resolve("pfsearch_parsed.json")
    def json = JsonOutput.toJson(profilesMatches)
    new File(outputFilePath.toString()).write(json)
//     pfsearch_parser.py \
//         ${pfsearch_out} \
//         ${pfsearch_out}-filtered.json \
//         ${blacklist_file}
}
