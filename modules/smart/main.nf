import java.nio.file.Files
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.StandardCopyOption
import java.util.zip.ZipEntry
import java.util.zip.ZipFile
import java.util.zip.ZipOutputStream
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.FastaFile
import uk.ac.ebi.interpro.HMMER2
import uk.ac.ebi.interpro.HMMER3

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
    def sequences = FastaFile.parse(fasta)  // [md5: sequence]
    def matches = HMMER3.parseOutput(hmmseach_out, "SMART")

    // Map model accessions to seq Ids -> [modelAcc: [seqIds]]
    def model2seqs = [:].withDefault { [] }
    matches.each { seqId, seqMatches ->
        seqMatches.each { modelAcc, _ ->
            model2seqs[modelAcc] << seqId
        }
    }

    def zipFile = task.workDir.resolve("sequences.zip")
    Files.newOutputStream(zipFile).withCloseable { os ->
        new ZipOutputStream(os).withCloseable { zip ->
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

                ZipEntry entry = new ZipEntry(fastaFile.fileName.toString())
                zip.putNextEntry(entry)
                Files.copy(fastaFile, zip)
                zip.closeEntry()
                Files.delete(fastaFile);
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

    unzip -d /tmp/fasta ${seqs_zip}
    touch /tmp/hmmpfam.out

    for fasta in /tmp/fasta/*.faa; do
        name=\${fasta##*/}   # strip directory
        name=\${name%.faa}   # strip extension
        hmm="${hmmdir}/\${name}.hmm"
        if [[ -f "\$hmm" ]]; then
            hmmpfam --acc -A 0 -E 0.01 -Z 350000 --cpu ${task.cpus} "\$hmm" "\$fasta" >> /tmp/hmmpfam.out
            zip -jq /tmp/smart.zip "\$hmm"
        fi
    done

    zip -jq /tmp/smart.zip /tmp/hmmpfam.out
    rm -r /tmp/fasta /tmp/hmmpfam.out
    mv /tmp/smart.zip smart.zip
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
    Files.newInputStream(fasta_zip).withCloseable { is ->
        Files.copy(is, tmp_zip, StandardCopyOption.REPLACE_EXISTING)
    }

    new ZipFile(tmp_zip.toFile()).withCloseable { zip ->
        zip.entries().each { entry ->
            if (entry.isDirectory()) {
                // Should never happen
                return
            }

            def fileName = Paths.get(entry.name).fileName.toString()
            Path extracted = task.workDir.resolve(fileName)

            zip.getInputStream(entry).withCloseable { input ->
                Files.copy(input, extracted, StandardCopyOption.REPLACE_EXISTING)
            }

            sequences = sequences + FastaFile.parse(extracted)
            Files.deleteIfExists(extracted)
        }
    }

    Files.newInputStream(search_zip).withCloseable { is ->
        Files.copy(is, tmp_zip, StandardCopyOption.REPLACE_EXISTING)
    }

    new ZipFile(tmp_zip.toFile()).withCloseable { zip ->
        zip.entries().each { entry ->
            if (entry.isDirectory()) {
                // Should never happen
                return
            }

            def fileName = Paths.get(entry.name).fileName.toString()
            Path extracted = task.workDir.resolve(fileName)

            zip.getInputStream(entry).withCloseable { input ->
                Files.copy(input, extracted, StandardCopyOption.REPLACE_EXISTING)
            }

            if (fileName == "hmmpfam.out") {
                hmmpfam_out = extracted
            } else {
                def (accession, length) = HMMER2.parseHMM(extracted)
                hmm_lengths[accession] = length
                Files.deleteIfExists(extracted)
            }
        }
    }

    Files.deleteIfExists(tmp_zip)
   
    def matches = HMMER2.parseOutput(hmmpfam_out, hmm_lengths, "SMART")

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
