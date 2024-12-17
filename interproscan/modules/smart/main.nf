import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process SEARCH_SMART {
    label 'medium', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmbindb

    output:
    tuple val(meta), path("hmmpfam.out")

    script:
    """
    /opt/hmmer2/bin/hmmpfam \
        --acc -A 0 -E 0.01 -Z 350000 \
        ${hmmbindb} ${fasta} > hmmpfam.out
    """
}

process PARSE_SMART {
    label 'small'

    input:
    tuple val(meta), val(hmmpfam_out), val(seq_json)
    val hmmtxtdb

    output:
    tuple val(meta), path("smart.json")

    exec:
    def jsonFile = new File(seq_json.toString())
    def jsonSlurper = new JsonSlurper()
    def sequences = jsonSlurper.parse(jsonFile)
        .collectEntries{ seqId, obj ->
            if (obj instanceof List) { // nucleotide sequences case
                obj.collectEntries { seq ->
                    [(seq.id): FastaSequence.fromMap(seq)]
                }
            } else {
                [(seqId): FastaSequence.fromMap(obj)]
            }
        }

    def hmmLengths = HMMER2.parseHMM(hmmtxtdb.toString())
    def matches = HMMER2.parseOutput(hmmpfam_out.toString(), hmmLengths)

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
