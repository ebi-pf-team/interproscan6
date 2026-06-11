process REPORT_NO_MATCHES {
    label    'mem_low'
    label    'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(fasta), val(jsons)

    output:
    tuple val(meta), path("no_matches.json")

    exec:
    def seqs_with_matches = [] as Set
    jsons.each { json ->
        def matches = json.newReader().withCloseable { reader ->
            new com.fasterxml.jackson.databind.ObjectMapper().readValue(reader, Map)
        }
        matches.each { seq_md5, _seq_matches ->
            seqs_with_matches.add(seq_md5)
        }
    }

    def sequences = uk.ac.ebi.interpro.FastaFile.parse(fasta)
    def no_matches = [:]
    sequences.each { seq_md5, _sequence ->
        if (!seqs_with_matches.contains(seq_md5)) {
            no_matches[seq_md5] = [:]
        }
    }
    def filepath = task.workDir.resolve("no_matches.json")
    filepath.text = groovy.json.JsonOutput.toJson(no_matches)
}
