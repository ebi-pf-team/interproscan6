import groovy.json.JsonOutput

process SEARCH_SMART {
    label 'medium', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmbindb

    output:
    tuple val(meta), path("hmmpfam.out"), path(fasta)

    script:
    """
    /opt/hmmer2/bin/hmmpfam \
        --acc -A 0 -E 0.01 -Z 350000 \
        ${hmmbindb} ${fasta} > hmmpfam.out
    """
}

process PARSE_SMART {
    label 'run_locally'

    input:
    tuple val(meta), val(hmmpfam_out), val(fasta)
    val hmmtxtdb

    output:
    tuple val(meta), path("smart.json")

    exec:
    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]

    def hmmLengths = HMMER2.parseHMM(hmmtxtdb.toString())
    def matches = HMMER2.parseOutput(hmmpfam_out.toString(), hmmLengths, "SMART")

    String tyrKinaseAccession = "SM00219"
    def tyrKinasePattern = ~/.*HRD[LIV][AR]\w\wN.*/
    String serThrKinaseAccession = "SM00220"
    def serThrKinasePattern = ~/.*D[LIVM]K\w\wN.*/

    matches.collectEntries { seqId, models ->
        def filteredModels = [:]

        if (models.containsKey(tyrKinaseAccession) && models.containsKey(serThrKinaseAccession)) {
            /*
                If both Tyrosine Kinase and Serine-Threonine Kinase
                have a hit against the same sequence,
                we need to perform an additional check before selecting them
            */
            String sequence = sequences[seqId].sequence
            boolean tyrKinaseOK = (sequence ==~ tyrKinasePattern)
            boolean serThrKinaseOK = (sequence ==~ serThrKinasePattern)

            models.each { modelAccession, match ->
                if (modelAccession != tyrKinaseAccession &&
                    modelAccession != serThrKinaseAccession) {
                    filteredModels[modelAccession] = match
                } else if (modelAccession == tyrKinaseAccession && tyrKinaseOK) {
                    filteredModels[modelAccession] = match
                } else if (modelAccession == serThrKinaseAccession && serThrKinaseOK) {
                    filteredModels[modelAccession] = match
                }
            }
        } else {
            filteredModels = models
        }

        [ (seqId): filteredModels ]
    }

    def outputFilePath = task.workDir.resolve("smart.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
