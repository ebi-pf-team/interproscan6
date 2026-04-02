import groovy.io.FileType
import groovy.json.JsonOutput
import groovy.json.JsonSlurper
import uk.ac.ebi.interpro.Location
import uk.ac.ebi.interpro.Match
import uk.ac.ebi.interpro.Signature
import uk.ac.ebi.interpro.SignatureLibraryRelease


process RUN_PSSCAN {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the
    pftools developers. It automates running pfscan for all provided patterns and
    includes post-processing of the hits.
    */
    label 'mem_min', 'time_medium', 'ips6_container'

    input:
        tuple val(meta), path(fasta)
        tuple path(dat), path(eval)

    output:
        tuple val(meta), path("ps_scan.out")

    script:
    """
    ps_scan.pl \
        ${fasta} \
        -d ${dat} \
        --pfscan pfscanV3 \
        -b ${eval} \
        -r -s -o ipro > ps_scan.out
    """
}


process PARSE_PSSCAN {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
        tuple val(meta), val(ps_scan_out)

    output:
        tuple val(meta), path("ps_scan.json")

    exec:
    Map<String, Map<String, Match>> patternsMatches = [:]
    SignatureLibraryRelease library = new SignatureLibraryRelease("PROSITE patterns", null)
    ps_scan_out.eachLine { line ->
        line = line.trim()
        if (!line || line.startsWith("pfscanV3 is not meant to be used with a single profile")) {
            return
        }
        List<String> matchInfo = line.split('\t')
        if (matchInfo.size() < 9) {
            return
        }
        String seqId = matchInfo[0]
        String modelAccession = matchInfo[2]
        int start = matchInfo[3].toInteger()
        int end = matchInfo[4].toInteger()

        List<String> matchDetails = matchInfo[8].split(';')
        String name = matchDetails[0].trim()
        String level = matchDetails[1].trim()
        if (!level.startsWith("LevelTag") || !level.contains("0")) {
            return // skipping non-strong matches
        } else {
            level = "STRONG"
        }
        String alignment = matchDetails[2].replaceAll('Sequence ', '').replaceAll('"', '').replaceAll('\\.', '').trim()
        String cigarAlignment = Match.encodeCigarAlignment(alignment)
        patternsMatches.computeIfAbsent(seqId) { [:] }
        Match matchObj = patternsMatches[seqId].computeIfAbsent(modelAccession) {
            new Match(modelAccession, new Signature(modelAccession, library))
        }
        Location location = new Location(start, end, level, alignment, cigarAlignment)
        matchObj.addLocation(location)
    }

    def filepath = task.workDir.resolve("ps_scan.json")
    filepath.text = JsonOutput.toJson(patternsMatches)
}

process RUN_PFSEARCH {
    label 'mem_min', 'time_medium', 'dynamic', 'ips6_container'

    input:
    tuple val(meta), path(fasta)
    path profiles_dir

    output:
    tuple val(meta), stdout

    script:
    """
    find ${profiles_dir}/ -type f | while read profile; do
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
    val blacklist_file

    output:
    tuple val(meta), path("pfsearch.json")

    exec:
    Map matches = [:]
    SignatureLibraryRelease library = new SignatureLibraryRelease(signature_library, null)
    def toSkip = []
    if (blacklist_file) {
        toSkip = blacklist_file.readLines()
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
    def filepath = task.workDir.resolve("pfsearch.json")
    filepath.text = JsonOutput.toJson(matches)
}

