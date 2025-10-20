import java.nio.file.Files
import groovy.json.JsonOutput

process PREFILTER_SMART {
    label 'small', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val hmmfile

    output:
    tuple val(meta), path("hmmsearch.out"), path(fasta)

    script:
    """
    hmmsearch \
        -E 100 --domE 100 --incE 100 --incdomE 100 --cpu ${task.cpus} \
        ${dirpath}/${hmmfile} ${fasta} > hmmsearch.out
    """
}

process PREPARE_SMART {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_out), val(fasta)
    val dirpath
    val hmmdir
    val chunk_size

    output:
    tuple val(meta), path(fastaFiles), val(smarts)

    exec:
    /* Only run seqs against a HMMER2 model where a HMMER3 match was found.
        Return tuple val(meta), val(fastaFiles), val(smarts) where the Nth item of
        fastaFiles corresponds to the sequences to scan against the Nth element of smarts. */
    def matches = HMMER3.parseOutput(hmmseach_out.toString(), "SMART")

    // Extract model accessions against which at least one sequence match was found
    smarts = matches.values()
        .collect { jsonMatches -> jsonMatches.keySet() }
        .flatten()
        .unique()
        .findAll { smartId ->
            new File("${dirpath.toString()}/${hmmdir}/${smartId}.hmm").exists()
        }

    // Map model accessions to seq Ids -> [modelAcc: [seqIds]]
    Map<String, List<String>> model2seqs = [:].withDefault { [] }
    matches.each { seqId, seqMatches ->
        seqMatches.each { modelAcc, _ ->
            model2seqs[modelAcc] << seqId
        }
    }

    // Build custom FASTA files with only seqs that matched the HMMER3 smart model
    Map<String, String> allSeqs = FastaFile.parse(fasta.toString())
    fastaFiles = []
    smarts.each { modelAcc ->
        seqIds = model2seqs[modelAcc]
        Path fastaPath = task.workDir.resolve("model_fasta_${modelAcc}.fasta")
        new File(fastaPath.toString()).withWriter('UTF-8') { writer ->
            seqIds.each { seqId ->
                seq = allSeqs[seqId]
                writer.writeLine(">${seqId}")
                for (int i = 0; i < seq.length(); i += 60) {
                    writer.writeLine(seq.substring(i, Math.min(i + 60, seq.length())))
                }
            }
        }
        fastaFiles.add( fastaPath )
    }
}

process SEARCH_SMART {
    label 'medium', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta), val(smarts)
    path dirpath
    val hmmdir

    output:
    tuple val(meta), path("hmmpfam.out"), path(fasta)

    script:
    def commands = ""
    if (fasta.size() == 0) {
        // Create an empty hmmpfam.out file if input is empty
        commands = "touch hmmpfam.out"
    } else {
        smarts.each { smartFile ->
            fasta.each { chunkFile ->
                String hmmFilePath = "${dirpath.toString()}/${hmmdir}/${smartFile}.hmm"  // reassign to a var so the cmd can run
                commands += "hmmpfam"
                commands += " --acc -A 0 -E 0.01 -Z 350000 --cpu ${task.cpus}"
                commands += " $hmmFilePath ${chunkFile} >> hmmpfam.out\n"
            }
        }
    }

    """
    ${commands}
    """
}

process PARSE_SMART {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmpfam_out), val(fasta)
    val dirpath
    val hmmdir

    output:
    tuple val(meta), path("smart.json")

    exec:
    // fasta may be a single file or multiple
    Map<String, String> sequences = [:] // [md5: sequence]
    if (fasta instanceof List) {
        fasta.each { fastaFile ->
            sequences = sequences + FastaFile.parse(fastaFile.toString())
        }
    } else if (Files.exists(fasta)) {
        sequences = FastaFile.parse(fasta.toString())
    }

    def hmmLengths = HMMER2.parseHMMs("${dirpath.toString()}/${hmmdir}")
    def matches = HMMER2.parseOutput(hmmpfam_out.toString(), hmmLengths, "SMART")

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
    def json = JsonOutput.toJson(filteredMatches)
    new File(outputFilePath.toString()).write(json)
}
