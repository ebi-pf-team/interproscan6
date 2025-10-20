import groovy.json.JsonOutput

import Match

process PREPROCESS_HAMAP {
    label 'mini', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val hmmfile

    output:
    tuple val(meta), path("hmmsearch.tab"), path(fasta)

    script:
    """
    hmmsearch \
        -E 100 --domE 100 --incE 100 --incdomE 100 \
        --cpu ${task.cpus} \
        --tblout hmmsearch.tab \
        ${dirpath}/${hmmfile} ${fasta} > /dev/null
    """
}

process PREPARE_HAMAP {
    label 'tiny'
    executor 'local'

    input:
    tuple val(meta), val(hmmsearch_tab), val(fasta)
    val dirpath
    val profiles_dir

    output:
    tuple val(meta), val(profiles), path(fasta_files)

    exec:
    Map<String, String> sequences = FastaFile.parse(fasta.toString())  // [md5: sequence]

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
    matches
        .each { query, targets ->
            Path prfPath = file("${dirpath.toString()}/${profiles_dir}/${query}.prf")
            File prfFile = new File(prfPath.toString())
            if (prfFile.exists()) {
                Path fastaPath = task.workDir.resolve("${query}.fa")
                new File(fastaPath.toString()).withWriter('UTF-8') { writer ->
                    targets.each { seqId ->
                        String seq = sequences[seqId]
                        writer.writeLine(">${seqId}")
                        writer.writeLine(fmtSequence(seq))
                    }
                }

                profiles.add( query )
                fasta_files.add( fastaPath )
            }
        }
}

process RUN_HAMAP {
    label 'mini', 'ips6_container'

    input:
    tuple val(meta), val(profiles), path(fasta_files)
    path dirpath
    val profiles_dir

    output:
    tuple val(meta), stdout

    script:
    def commands = ""

    [profiles, fasta_files]
        .transpose()
        .each { profile, fasta ->
            def profilePath = "${dirpath.toString()}/${profiles_dir}/${profile}.prf"
            commands += "pfsearchV3 -f -o 7 ${profilePath} ${fasta}\n"
        }

    """
    ${commands}
    """
}

process PARSE_HAMAP {
    label 'tiny'
    executor 'local'

    input:
    tuple val(meta), val(pfsearch_out)

    output:
    tuple val(meta), path("hamap.json")

    exec:
    SignatureLibraryRelease library = new SignatureLibraryRelease("HAMAP", null)
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
        String cigarAlignment = Match.encodeCigarAlignment(alignment)

        if (matches.containsKey(sequenceId) && matches[sequenceId].containsKey(modelAccession)) {
            match = matches[sequenceId][modelAccession]
        } else {
            match = new Match(modelAccession, new Signature(modelAccession, library))
            matches.computeIfAbsent(sequenceId, { [:] })
            matches[sequenceId][modelAccession] = match
        }
        Location location = new Location(start, end, score, alignment, cigarAlignment)
        match.addLocation(location)
    }

    def outputFilePath = task.workDir.resolve("hamap.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}

def fmtSequence(String sequence) {
    /* Use a stringBuild for efficiency, this stops a new str being created
    with each addition of a new line char.*/
    StringBuilder sb = new StringBuilder();
    for (int i = 0; i < sequence.length(); i += 60) {
        int j = Math.min(i + 60, sequence.length());
        sb.append(sequence, i, j);
        if (j < sequence.length()) {
            sb.append('\n');
        }
    }
    return sb.toString()
}
