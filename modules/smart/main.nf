process PREFILTER_SMART {
    label     'mem_min', 'time_veryshort', 'dynamic'
    container 'interpro/smart:1.0'

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

    output:
    tuple val(meta), path("sequences.zip")

    exec:
    def sequences = uk.ac.ebi.interpro.FastaFile.parse(fasta)  // [md5: sequence]
    def matches = uk.ac.ebi.interpro.HMMER3.parseOutput(hmmseach_out, "SMART")

    // Map model accessions to seq Ids -> [modelAcc: [seqIds]]
    def model2seqs = [:].withDefault { [] }
    matches.each { seqId, seqMatches ->
        seqMatches.each { modelAcc, match ->
            model2seqs[modelAcc] << seqId
        }
    }

    def zipFile = task.workDir.resolve("sequences.zip")
    java.nio.file.Files.newOutputStream(zipFile).withCloseable { os ->
        new java.util.zip.ZipOutputStream(os).withCloseable { zip ->
            // Create a FASTA file for each profile to search with
            model2seqs.sort().each { modelAcc, seqIds ->
                def fastaFile = task.workDir.resolve("${modelAcc}.faa")
                fastaFile.withWriter { writer ->
                    seqIds.sort().each { seqId ->
                        def seq = sequences[seqId]
                        writer.writeLine(">${seqId}")
                        seq.eachMatch(/.{1,60}/) { writer.writeLine(it) }
                    }
                }

                def entry = new java.util.zip.ZipEntry(fastaFile.fileName.toString())
                zip.putNextEntry(entry)
                java.nio.file.Files.copy(fastaFile, zip)
                zip.closeEntry()
                java.nio.file.Files.delete(fastaFile);
            }
        }
    }
}

process SEARCH_SMART {
    label     'mem_min', 'time_veryshort', 'dynamic'
    container 'interpro/smart:1.0'

    input:
    tuple val(meta), path(seqs_zip)
    path hmmdir

    output:
    tuple val(meta), path(seqs_zip), path("smart.zip")

    script:
    """
    shopt -s nullglob

    tmpdir=\$(mktemp -d -p .)

    unzip -q ${seqs_zip} -d \${tmpdir}
    touch \${tmpdir}/hmmpfam.out

    for fasta in \${tmpdir}/*.faa; do
        name=\${fasta##*/}   # strip directory
        name=\${name%.faa}   # strip extension
        hmm="${hmmdir}/\${name}.hmm"
        if [[ -f "\$hmm" ]]; then
            hmmpfam --acc -A 0 -E 0.01 -Z 350000 --cpu ${task.cpus} "\$hmm" "\$fasta" >> \${tmpdir}/hmmpfam.out
            zip -jq \${tmpdir}/smart.zip "\$hmm"
        fi
    done

    zip -jq \${tmpdir}/smart.zip \${tmpdir}/hmmpfam.out
    mv \${tmpdir}/smart.zip smart.zip
    rm -r \${tmpdir}
    """
}

process PARSE_SMART {
    label    'mem_low', 'time_short'
    executor 'local'

    input:
    tuple val(meta), val(fasta_zip), val(search_zip)

    output:
    tuple val(meta), path("smart.json")

    exec:
    def sequences   = [:]
    def hmm_lengths = [:]
    def hmmpfam_out = null

    def tmp_zip = task.workDir.resolve("temp.zip")
    java.nio.file.Files.newInputStream(fasta_zip).withCloseable { is ->
        java.nio.file.Files.copy(is, tmp_zip, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }

    new java.util.zip.ZipFile(tmp_zip.toFile()).withCloseable { zip ->
        zip.entries().each { entry ->
            if (entry.isDirectory()) {
                // Should never happen
                return
            }

            def fileName = java.nio.file.Paths.get(entry.name).fileName.toString()
            def extracted = task.workDir.resolve(fileName)

            zip.getInputStream(entry).withCloseable { input ->
                java.nio.file.Files.copy(input, extracted, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            }

            sequences = sequences + uk.ac.ebi.interpro.FastaFile.parse(extracted)
            java.nio.file.Files.deleteIfExists(extracted)
        }
    }

    java.nio.file.Files.newInputStream(search_zip).withCloseable { is ->
        java.nio.file.Files.copy(is, tmp_zip, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
    }

    new java.util.zip.ZipFile(tmp_zip.toFile()).withCloseable { zip ->
        zip.entries().each { entry ->
            if (entry.isDirectory()) {
                // Should never happen
                return
            }

            def fileName = java.nio.file.Paths.get(entry.name).fileName.toString()
            def extracted = task.workDir.resolve(fileName)

            zip.getInputStream(entry).withCloseable { input ->
                java.nio.file.Files.copy(input, extracted, java.nio.file.StandardCopyOption.REPLACE_EXISTING)
            }

            if (fileName == "hmmpfam.out") {
                hmmpfam_out = extracted
            } else {
                def (accession, length) = uk.ac.ebi.interpro.HMMER2.parseHMM(extracted)
                hmm_lengths[accession] = length
                java.nio.file.Files.deleteIfExists(extracted)
            }
        }
    }

    java.nio.file.Files.deleteIfExists(tmp_zip)
   
    def matches = uk.ac.ebi.interpro.HMMER2.parseOutput(hmmpfam_out, hmm_lengths, "SMART")

    def tyrKinaseAccession = "SM00219"
    def tyrKinasePattern = ~/.*HRD[LIV][AR]\w\wN.*/
    def serThrKinaseAccession = "SM00220"
    def serThrKinasePattern = ~/.*D[LIVM]K\w\wN.*/

    def filteredMatches = matches.collectEntries { seqId, models ->
        def filteredModels = [:]

        if (models.containsKey(tyrKinaseAccession) && models.containsKey(serThrKinaseAccession)) {
            /*
                If both Tyrosine Kinase and Serine-Threonine Kinase
                have a hit against the same sequence,
                we need to perform an additional check before selecting them
            */
            def sequence = sequences[seqId]
            def tyrKinaseOK = (sequence ==~ tyrKinasePattern)
            def serThrKinaseOK = (sequence ==~ serThrKinasePattern)

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
    outputFilePath.text  = groovy.json.JsonOutput.toJson(filteredMatches)
}
