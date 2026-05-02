process SEARCH_SUPERFAMILY {
    label     'mem_min', 'time_short', 'dynamic'
    container 'interpro/hmmer:3.3'

    input:
    tuple val(meta), val(meta2), path(fasta)
    tuple path(dir), val(hmm), val(selfhits), val(cla), val(model_tab), val(pdbj95d)

    output:
    tuple val(meta), val(meta2), path("superfamily.out")

    script:
    """
    hmmscan \
        -E 10 -Z 15438 \
        --cpu ${task.cpus} \
        ${dir}/${hmm} ${fasta} > hmmscan.out

    ass3_single_threaded.pl \
        -e 0.0001 -t n -f 1 \
        -s ${dir}/${selfhits} \
        -r ${dir}/${cla} \
        -m ${dir}/${model_tab} \
        -p ${dir}/${pdbj95d} \
        ${fasta} \
        hmmscan.out \
        superfamily.out
    """
}

process PARSE_SUPERFAMILY {
    label    'mem_low', 'time_veryshort'
    executor 'local'

    input:
    tuple val(meta), val(meta2), val(superfamily_out)
    tuple val(dirpath), val(hmm), val(model_tab)

    output:
    tuple val(meta), val(meta2), path("superfamily.json")

    exec:
    def library = new uk.ac.ebi.interpro.SignatureLibraryRelease("SUPERFAMILY", null)
    def model2sf = [:]
    dirpath.resolve(model_tab).eachLine { line ->
        def fields = line.trim().split(/\t/)
        def modelId = fields[0]
        def superfamilyAccession = fields[1]
        assert !model2sf.containsKey(modelId)
        model2sf[modelId] = "SSF${superfamilyAccession}"
    }

    def model2length = [:]
    def modelAc = null
    def length = null
    dirpath.resolve(hmm).eachLine { line ->
        line = line.trim()
        if (line.startsWith('//')) {
            assert modelAc != null && length != null
            model2length[modelAc] = length
            modelAc = null
            length = null
        } else if (line.startsWith('N') && !modelAc) {
            def match = (line =~ ~/^NAME\s+(.+)$/)
            if (match) modelAc = match[0][1]
        } else if (line.startsWith('L') && !length) {
            def match = (line =~ ~/^LENG\s+([0-9]+)$/)
            if (match) length = match[0][1].toInteger()
        }
    }

    def matches = [:].withDefault { [:] }
    superfamily_out.eachLine { line ->
        line = line.trim()
        if (line) {
            def fields = line.split(/\s+/)
            assert fields.size() == 9
            def seqId = fields[0]
            def modelId = fields[1]
            if (modelId != "-") {
                def superfamilyAccession = model2sf[modelId]
                assert superfamilyAccession != null

                def regionsAsString = fields[2]
                def evalue = Double.parseDouble(fields[3])

                def regions = []
                regionsAsString.split(",").each { region ->
                    def boundaries = region.split("-")
                    assert boundaries.size() == 2
                    def start = boundaries[0].toInteger()
                    def end = boundaries[1].toInteger()
                    regions.add([start, end])
                }

                assert regions.size() >= 1

                // Sort by start/end
                regions = regions.sort { a, b ->
                    a[0] <=> b[0] ?: a[1] <=> b[1]
                }

                def start = regions[0][0]
                def end = regions.collect { r -> r[1] }.max()
                def hmmLength = model2length[modelId]
                def fragments = []
                if (regions.size() > 1) {
                    regions.eachWithIndex { obj, idx ->
                        def (fragStart, fragEnd) = obj
                        def dcStatus
                        if (idx == 0) {
                            dcStatus = "C_TERMINAL_DISC"
                        } else if (idx == regions.size() - 1) {
                            dcStatus = "N_TERMINAL_DISC"
                        } else {
                            dcStatus = "NC_TERMINAL_DISC"
                        }
                        fragments.add(new uk.ac.ebi.interpro.LocationFragment(fragStart, fragEnd, dcStatus))
                    }

                } else {
                    def (fragStart, fragEnd) = regions[0]
                    fragments.add(new uk.ac.ebi.interpro.LocationFragment(fragStart, fragEnd, "CONTINUOUS"))
                }

                def location = new uk.ac.ebi.interpro.Location(start, end, hmmLength, evalue, fragments)
                def match = matches[seqId][modelId]
                if (match == null) {
                    match = new uk.ac.ebi.interpro.Match(modelId, new uk.ac.ebi.interpro.Signature(superfamilyAccession, library))
                    match.addLocation(location)
                    matches[seqId][modelId] = match
                } else {
                    match.addLocation(location)
                }
            }
        }
    }

    def filepath = task.workDir.resolve("superfamily.json")
    filepath.text  = groovy.json.JsonOutput.toJson(matches)
}
