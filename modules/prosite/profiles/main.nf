import groovy.io.FileType
import groovy.json.JsonOutput
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease

process RUN_PFSEARCH {
    label 'mem_min', 'time_medium', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path dirpath
    val profiles_dir

    output:
    tuple val(meta), stdout

    script:
    """
    find ${dirpath}/${profiles_dir} -type f | while read profile; do
        pfsearchV3 -f -o 7 -t ${task.cpus} "\${profile}" "${fasta}"
    done
    """
}

process PARSE_PFSEARCH {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(pfsearch_out)
    val signature_library
    val dirpath
    val blacklist_file

    output:
    tuple val(meta), path("pfsearch.json")

    exec:
    Map matches = [:]
    SignatureLibraryRelease library = new SignatureLibraryRelease(signature_library, null)
    def toSkip = []
    if (dirpath && blacklist_file) {
        toSkip = new File("${dirpath.toString()}/${blacklist_file}").readLines()
    }

    pfsearch_out.eachLine { line ->
        def fields = line.split()
        assert fields.size() == 10
        String modelAccession = fields[0].split("\\|")[0]
        if (toSkip && (modelAccession in toSkip)) {
            return // skip flagged accessions
        }

        String seqId = fields[3]
        int start = fields[4].toInteger()
        int end = fields[5].toInteger()
        Double score = Double.parseDouble(fields[7])
        String alignment = fields[9]
        String cigarAlignment = Match.encodeCigarAlignment(alignment)

        matches.computeIfAbsent(seqId) { [:] }
        Match matchObj = matches[seqId].computeIfAbsent(modelAccession) {
            new Match(modelAccession, new Signature(modelAccession, library))
        }
        Location location = new Location(start, end, score, alignment, cigarAlignment)
        matchObj.addLocation(location)
    }
    def outputFilePath = task.workDir.resolve("pfsearch.json")
    def json = JsonOutput.toJson(matches)
    new File(outputFilePath.toString()).write(json)
}
