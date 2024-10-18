import groovy.json.JsonOutput
import groovy.json.JsonSlurper

process PREPROCESS_HAMAP {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb

    output:
    tuple val(meta), path("hmmsearch.tab")

    script:
    """
    /opt/hmmer3/bin/hmmsearch \
        -E 100 --domE 100 --incE 100 --incdomE 100 \
        --cpu ${task.cpus} \
        --tblout hmmsearch.tab \
        ${hmmdb} ${fasta} > /dev/null
    """
}

process PREPARE_HAMAP {
    label 'analysis_parser'

    input:
    tuple val(meta), val(hmmsearch_tab)
    val seq_json
    val profile_dir

    output:
    tuple val(meta), val(profiles), path(fasta_files)

    exec:
    def jsonFile = new File(seq_json.toString())
    def jsonSlurper = new JsonSlurper()
    def sequences = jsonSlurper.parse(jsonFile)
        .collectEntries{ seqId, obj ->
            [ seqId, FastaSequence.fromMap(obj) ]
        }

    // Find profiles with matches
    def matches = [:]
    file(hmmsearch_tab.toString()).eachLine { line ->
        if (line[0] != "#") {
            def fields = line.trim().split(/\s+/, 19)
            String target = fields[0]
            String query = fields[2]

            if (matches.containsKey(query)) {
                matches[query].add(target)
            } else {
                matches[query] = [target]
            }
        }
    }

    // Do not use `def` so lists are globally scoped and can be used in the `output` block
    profiles = []
    fasta_files = []

    // Create a FASTA file for each profile to search with
    matches = matches
        .each { query, targets ->
            Path prfPath = file("${profile_dir.toString()}/${query}.prf")
            File prfFile = new File(prfPath.toString())
            if (prfFile.exists()) {
                Path fastaPath = task.workDir.resolve("${query}.fa")
                new File(fastaPath.toString()).withWriter('UTF-8') { writer ->
                    targets.each { seqId ->
                        FastaSequence seq = sequences[seqId]
                        writer.writeLine(">${seqId}")
                        writer.writeLine(seq.sequence)
                    }
                }

                profiles.add( query )
                fasta_files.add( fastaPath )
            }
        }
}

process RUN_HAMAP {
    label 'analysis_parser'

    input:
    tuple val(meta), val(profiles), path(fasta_files)
    path profile_dir

    output:
    tuple val(meta), stdout

    script:
    def commands = ""

    [profiles, fasta_files]
        .transpose()
        .each { profile, fasta ->
            def profilePath = "${profile_dir.toString()}/${profile}.prf"
            commands += "/opt/pftools/pfsearchV3 -f -o 7 ${profilePath} ${fasta}\n"
        }

    """
    ${commands}
    """
}

process PARSE_HAMAP {
    label 'analysis_parser'

    input:
    tuple val(meta), val(pfsearch_out)

    output:
    tuple val(meta), path("hamap.json")

    exec:
    def matches = [:]
    pfsearch_out.eachLine { line ->
        def fields = line.split()
        assert fields.size() == 10
        String modelAccession = fields[0].split("\\|")[0]
        String sequenceId = fields[3]
        int start = fields[4].toInteger()
        int end = fields[5].toInteger()
        Double score = Double.parseDouble(fields[7])
        String alignment = fields[9]

        if (matches.containsKey(sequenceId)) {
            match = matches[sequenceId]
        } else {
            match = new Match(modelAccession)
            matches[sequenceId] = [:]
            matches[sequenceId][modelAccession] = match
        }

        Location location = new Location(start, end, score, alignment)
        match.addLocation(location)
    }

    def outputFilePath = task.workDir.resolve("hamap.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
