import groovy.json.JsonSlurper
import groovy.json.JsonOutput

process PREPARE_INTERPRO_N {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(fasta)

    output:
    tuple val(meta), path("sequences.tsv")

    exec:
    def MAX_LENGTH = 2047
    def outputFilePath = task.workDir.resolve("sequences.tsv")
    new File(outputFilePath.toString()).withWriter("UTF-8") { writer ->
        writer.writeLine("accession\tsequence")

        def identifier = ""
        def sequence = new StringBuilder()

        new File(fasta.toString()).eachLine { line ->
            line = line.trim()
            if (line.startsWith(">")) {
                if (identifier && sequence.length() <= MAX_LENGTH) {
                    writer.writeLine("${identifier}\t${sequence.toString()}")
                }
                identifier = line.substring(1).split("\\s+")[0]
                sequence.setLength(0)
            } else {
                sequence.append(line)
            }
        }

        if (identifier && sequence.length() <= MAX_LENGTH) {
            writer.writeLine("${identifier}\t${sequence.toString()}")
        }
    }
}

process RUN_INTERPRO_N_CPU {
    label 'large', 'interpro_n_container'

    input:
    tuple val(meta), path(tsv)
    path ck_dir
    val batch_size

    output:
    tuple val(meta), path("*json")

    script:
    """
    WORKDIR=\$(pwd)
    INPUT=\$(readlink -e ${tsv})
    CHECKPOINT=\$(readlink -e ${ck_dir})
    cd /workdir
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --n_seq_per_output_file 5 \
        --output_dir \${WORKDIR}
    """
}

process RUN_INTERPRO_N_GPU {
    label 'large', 'interpro_n_container', 'use_gpu'

    input:
    tuple val(meta), path(tsv)
    path ck_dir
    val batch_size

    output:
    tuple val(meta), path("*json")

    script:
    """
    WORKDIR=\$(pwd)
    INPUT=\$(readlink -e ${tsv})
    CHECKPOINT=\$(readlink -e ${ck_dir})
    cd /workdir
    python3 -m interpro_n.main \
        --input_dir \${INPUT} \
        --checkpoint \${CHECKPOINT} \
        --batch_size ${batch_size} \
        --n_seq_per_output_file 5 \
        --output_dir \${WORKDIR}
    """
}

process PARSE_INTERPRO_N {
    label    'tiny'
    executor 'local'

    input:
    tuple val(meta), val(json_files)
    val applications

    output:
    tuple val(meta), path("interpro-n.json")

    exec:
    def results = [:]
    def jsonSlurper = new JsonSlurper()

    json_files.each { jsonFile ->
        def file = new File(jsonFile.toString())
        def output = jsonSlurper.parse(file)
        output.each { item ->
            def seqId = item.accession
            def matches = [:]

            item.match.each { rawMatch ->
                def sigLib = rawMatch.signature.signatureLibraryRelease
                def sigLibName = sigLib.library.toLowerCase().replaceAll("[- ]", "")
                if (!applications.contains(sigLibName)) {
                    return
                }

                assert rawMatch.locations.size() == 1
                def loc = rawMatch.locations[0]
                def accession = rawMatch.signature.accession
                // use a namespaced key to avoid clashes with traditional applications
                def matchKey = "InterPro-N::${accession}".toString()
                
                results.computeIfAbsent(seqId) { [:] }
                results[seqId].computeIfAbsent(matchKey) {   
                    SignatureLibraryRelease library = new SignatureLibraryRelease(sigLib.library, sigLib.version)
                    Signature signature = new Signature(accession, library)
                    Match match = new Match(accession)
                    match.signature = signature
                    match.source = "InterPro-N"
                    match
                }

                // collect and sort raw fragments
                def rawFragments = loc["location-fragments"].collect { f ->
                    [f.start, f.end]
                }.sort { a, b -> a[0] <=> b[0] }

                // transform to LocationFragment objects
                List<LocationFragment> fragments = []
                if (rawFragments.size() > 1) {
                    rawFragments.eachWithIndex { f, i ->
                        String status
                        if (i == 0) {
                            status = "C_TERMINAL_DISC"
                        } else if (i + 1 < rawFragments.size()) {
                            status = "NC_TERMINAL_DISC"
                        } else {
                            status = "N_TERMINAL_DISC"
                        }
                        fragments << new LocationFragment(f[0], f[1], status)
                    }
                } else {
                    def f = rawFragments[0]
                    fragments << new LocationFragment(f[0], f[1], "CONTINUOUS")
                }

                results[seqId][matchKey].addLocation(
                    new Location(loc.start, loc.end, loc.score as double, fragments)
                )
            }
        }
    }

    def json = JsonOutput.toJson(results)
    def outputFilePath = task.workDir.resolve("interpro-n.json")
    new File(outputFilePath.toString()).write(json)
}
