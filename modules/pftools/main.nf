process RUN_PSSCAN {
    /*
    The ps_scan.pl script is a wrapper for the pfscan tool that is provided by the
    pftools developers. It automates running pfscan for all provided patterns and
    includes post-processing of the hits.
    */
    label     'mem_min', 'time_medium'
    container 'interpro/pftools:3.2.12'

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
    def patternsMatches = [:]
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("PROSITE patterns", null)
    ps_scan_out.eachLine { line ->
        line = line.trim()
        if (!line || line.startsWith("pfscanV3 is not meant to be used with a single profile")) {
            return
        }
        def matchInfo = line.split('\t')
        if (matchInfo.size() < 9) {
            return
        }
        def seqId = matchInfo[0]
        def modelAccession = matchInfo[2]
        def start = matchInfo[3].toInteger()
        def end = matchInfo[4].toInteger()

        def matchDetails = matchInfo[8].split(';')
        def name = matchDetails[0].trim()
        def level = matchDetails[1].trim()
        if (!level.startsWith("LevelTag") || !level.contains("0")) {
            return // skipping non-strong matches
        } else {
            level = "STRONG"
        }
        def alignment = matchDetails[2].replaceAll('Sequence ', '').replaceAll('"', '').replaceAll('\\.', '').trim()
        def cigarAlignment = uk.ac.ebi.interpro.Match.encodeCigarAlignment(alignment)
        patternsMatches.computeIfAbsent(seqId) { [:] }
        def matchObj = patternsMatches[seqId].computeIfAbsent(modelAccession) {
            new uk.ac.ebi.interpro.Match(modelAccession, new uk.ac.ebi.interpro.Signature(modelAccession, library))
        }
        def location = new uk.ac.ebi.interpro.Location(start, end, level, alignment, cigarAlignment)
        matchObj.addLocation(location)
    }

    def filepath = task.workDir.resolve("ps_scan.json")
    filepath.text  = groovy.json.JsonOutput.toJson(patternsMatches)
}

process RUN_PFSEARCH {
    label     'mem_min', 'time_medium', 'dynamic'
    container 'interpro/pftools:3.2.12'

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
    def matches = [:]
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease(signature_library, null)
    def toSkip = []
    if (blacklist_file) {
        toSkip = blacklist_file.readLines()
    }

    pfsearch_out.eachLine { line ->
        def fields = line.split()
        assert fields.size() == 10
        def modelAccession = fields[0].split("\\|")[0]
        if (toSkip && (modelAccession in toSkip)) {
            return // skip flagged accessions
        }

        def seqId = fields[3]
        def start = fields[4].toInteger()
        def end = fields[5].toInteger()
        def score = Double.parseDouble(fields[7])
        def alignment = fields[9]
        def cigarAlignment = uk.ac.ebi.interpro.Match.encodeCigarAlignment(alignment)

        matches.computeIfAbsent(seqId) { [:] }
        def matchObj = matches[seqId].computeIfAbsent(modelAccession) {
            new uk.ac.ebi.interpro.Match(modelAccession, new uk.ac.ebi.interpro.Signature(modelAccession, library))
        }
        def location = new uk.ac.ebi.interpro.Location(start, end, score, alignment, cigarAlignment)
        matchObj.addLocation(location)
    }
    def filepath = task.workDir.resolve("pfsearch.json")
    filepath.text  = groovy.json.JsonOutput.toJson(matches)
}

