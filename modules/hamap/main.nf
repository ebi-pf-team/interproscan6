import groovy.json.JsonOutput
import uk.ac.ebi.interpro.Location  
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease

process RUN_HAMAP {
    label 'mem_min', 'time_veryshort', 'dynamic', 'ips6_container'

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

process PARSE_HAMAP {
    label    'mem_low', 'time_veryshort'
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
