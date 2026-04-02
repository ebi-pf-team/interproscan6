import java.nio.file.Files
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.HMMER2
import uk.ac.ebi.interpro.HMMER3

process PREFILTER_SMART {
    label 'mem_min', 'time_veryshort', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path hmmfile

    output:
    tuple val(meta), path("hmmsearch.out"), path(fasta)

    script:
    """
    hmmsearch \
        -E 100 --domE 100 --incE 100 --incdomE 100 --cpu ${task.cpus} \
        ${hmmfile} ${fasta} > hmmsearch.out
    """
}

process PREPARE_SMART {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_out), val(fasta)
    val hmmdir

    output:
    tuple val(meta), path(fastas), val(models)

    exec:
    /* Only run seqs against a HMMER2 model where a HMMER3 match was found.
        Return tuple val(meta), val(fastaFiles), val(smarts) where the Nth item of
        fastaFiles corresponds to the sequences to scan against the Nth element of smarts. */
    Map<String, String> sequences = FastaFile.parse(fasta)  // [md5: sequence]
    Map<String, Map> matches = HMMER3.parseOutput(hmmseach_out, "SMART")

    // Map model accessions to seq Ids -> [modelAcc: [seqIds]]
    Map<String, List<String>> model2seqs = [:].withDefault { [] }
    matches.each { seqId, seqMatches ->
        seqMatches.each { modelAcc, _ ->
            model2seqs[modelAcc] << seqId
        }
    }

    fastas = []
    models = []
    // Create a FASTA file for each profile to search with
    model2seqs.each { modelAcc, seqIds ->
        def mdlFile = hmmdir.resolve("${modelAcc}.hmm")
        if (!mdlFile.exists()) return
        def fasta = task.workDir.resolve("${modelAcc}.fa")
        fasta.withWriter('UTF-8') { writer ->
            seqIds.each { seqId ->
                String seq = sequences[seqId]
                writer.writeLine(">${seqId}")
                seq.eachMatch(/.{1,60}/) { writer.writeLine(it) }
            }
        }
        fastas << fasta
        models << modelAcc
    }
}

process SEARCH_SMART {
    label 'mem_min', 'time_veryshort', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fastas), val(models)
    path hmmdir

    output:
    tuple val(meta), path(fastas), path("hmmpfam.out")

    script:
    def commands = ""
    [fastas, models]
        .transpose()
        .each { entry -> 
            def fasta = entry[0]
            def model = entry[1]
            def hmm = hmmdir.resolve("${model}.hmm")
            commands += "hmmpfam"
            commands += " --acc -A 0 -E 0.01 -Z 350000 --cpu ${task.cpus}"
            commands += " ${hmm} ${fasta} >> hmmpfam.out\n"
        }

    """
    touch hmmpfam.out
    ${commands}
    """
}

process PARSE_SMART {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fastas), val(hmmpfam_out)
    val hmmdir

    output:
    tuple val(meta), path("smart.json")

    exec:
    // fasta may be a single file or multiple
    def sequences = [:] // [md5: sequence]
    fastas.each { fasta ->
        sequences = sequences + FastaFile.parse(fasta)
    }

    def hmmLengths = HMMER2.parseHMMs(hmmdir)
    def matches = HMMER2.parseOutput(hmmpfam_out, hmmLengths, "SMART")

    String tyrKinaseAccession = "SM00219"
    def tyrKinasePattern = ~/.*HRD[LIV][AR]\w\wN.*/
    String serThrKinaseAccession = "SM00220"
    def serThrKinasePattern = ~/.*D[LIVM]K\w\wN.*/

    def filteredMatches = matches.collectEntries { seqId, models ->
        def filteredModels = [:]

        if (models.containsKey(tyrKinaseAccession) && models.containsKey(serThrKinaseAccession)) {
            /*
                If both Tyrosine Kinase and Serine-Threonine Kinase
                have a hit against the same sequence,
                we need to perform an additional check before selecting them
            */
            String sequence = sequences[seqId]
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
    outputFilePath.text = JsonOutput.toJson(filteredMatches)
}
