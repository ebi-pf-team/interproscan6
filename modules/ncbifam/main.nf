import groovy.json.JsonOutput
import uk.ac.ebi.interpro.HMMER3

process PARSE_NCBIFAM {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(hmmseach_out)

    output:
    tuple val(meta), path("ncbifam.json")

    exec:
    def hmmer_matches = HMMER3.parseOutput(hmmseach_out, "NCBIFAM")
    hmmer_matches = hmmer_matches.collectEntries { seqId, matches ->
        [seqId, matches.collectEntries { modelAccession, match ->
            def updatedModelAccession = modelAccession.split("\\.")[0]
            match.modelAccession = updatedModelAccession
            match.signature.accession = updatedModelAccession
            [(updatedModelAccession): match]
        }]
    }

    def filepath = task.workDir.resolve("ncbifam.json")
    filepath.text = JsonOutput.toJson(hmmer_matches)
}
