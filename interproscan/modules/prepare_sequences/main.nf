import groovy.json.JsonOutput

process PREPARE_NUCLEIC_SEQUENCES {
    label 'io'

    input:
    // n_fasta:     original nucleic sequences file
    // p_fasta:     translated protein sequences file
    tuple val(n_fasta), val(p_fasta)

    output:
    val "sequences.json"

    exec:
    def sequences = FastaFile.parse(n_fasta.toString())
    def ntSequences = [:]
    sequences.each { seq ->
        ntSequences[seq.id] = seq
    }

    sequences = FastaFile.parse(p_fasta.toString())
    def output = [:]
    sequences.each { seq ->
        def match = seq.description =~ /source=(\S+)/
        assert match
        def source = match[0][1]
        def ntSeq = ntSequences[source]
        seq.translatedFrom = ntSeq
        output[seq.md5] = seq
    }
    def outputPath = task.workDir.resolve("sequences.json")
    def json = JsonOutput.toJson(output)
    new File(outputPath.toString()).write(json)
}


process PREPARE_PROTEIN_SEQUENCES {
    label 'io'

    input:
    val fasta

    output:
    val "sequences.json"

    exec:
    def output = [:]
    def sequences = FastaFile.parse(fasta.toString())
    sequences.each { seq ->
        output[seq.md5] = seq
    }
    def outputPath = task.workDir.resolve("sequences.json")
    def json = JsonOutput.toJson(output)
    new File(outputPath.toString()).write(json)
}
