import groovy.json.JsonOutput

process SEARCH_SUPERFAMILY {
    label 'hmmer_runner'

    input:
    tuple val(meta), path(fasta)
    path hmmdb
    path selfhits
    path cla
    path model
    path pdbj95d

    output:
    tuple val(meta), path("superfamily.out")

    script:
    """
    /opt/hmmer3/bin/hmmpress ${hmmdb}

    /opt/hmmer3/bin/hmmscan \
        -E 10 -Z 15438 \
        --cpu ${task.cpus} \
        ${hmmdb} ${fasta} > hmmscan.out

    perl ${projectDir}/bin/superfamily/ass3_single_threaded.pl \
        -e 0.0001 -t n -f 1 \
        -s ${selfhits} \
        -r ${cla} \
        -m ${model} \
        -p ${pdbj95d} \
        ${fasta} \
        hmmscan.out \
        superfamily.out
    """
}

process PARSE_SUPERFAMILY {
    input:
    tuple val(meta), val(superfamily_out)
    val model_tsv
    val hmmdb

    output:
    tuple val(meta), path("superfamily.json")

    exec:
    def model2sf = [:]
    file(model_tsv.toString()).eachLine { line ->
        def fields = line.trim().split(/\t/)
        String modelId = fields[0]
        String superfamilyAccession = fields[1]
        assert !model2sf.containsKey(modelId)
        model2sf[modelId] = "SSF${superfamilyAccession}"
    }

    def model2length = [:]
    String modelAc = null
    Integer length = null
    new File(hmmdb).eachLine { line ->
        line = line.trim()
        if (line.startsWith('//')) {
            assert modelAc != null && length != null
            model2length[modelAc] = length
            modelAc = length = null
        } else if (line.startsWith('N') && !modelAc) {
            def match = (line =~ ~/^NAME\s+(.+)$/)
            if (match) modelAc = match[0][1]
        } else if (line.startsWith('L') && !length) {
            def match = (line =~ ~/^LENG\s+([0-9]+)$/)
            if (match) length = match[0][1].toInteger()
        }
    }

    def matches = [:].withDefault { [:] }
    file(superfamily_out.toString()).eachLine { line ->
        line = line.trim()
        if (line) {
            def fields = line.split(/\s+/)
            assert fields.size() == 9
            String seqId = fields[0]
            String modelId = fields[1]
            if (modelId != "-") {
                String superfamilyAccession = model2sf[modelId]
                assert superfamilyAccession != null

                String regionsAsString = fields[2]
                Double evalue = Double.parseDouble(fields[3])

                def regions = []
                regionsAsString.split(",").each { region ->
                    def boundaries = region.split("-")
                    assert boundaries.size() == 2
                    int start = boundaries[0].toInteger()
                    int end = boundaries[1].toInteger()
                    regions.add([start, end])
                }

                assert regions.size() >= 1

                // Sort by start/end
                regions = regions.sort { a, b ->
                    a[0] <=> b[0] ?: a[1] <=> b[1]
                }

                int start = regions[0][0]
                int end = regions.collect { it[1] }.max()
                Integer hmmLength = model2length[modelId]
                List<LocationFragment> fragments = []
                if (regions.size() > 1) {
                    regions.eachWithIndex { obj, idx ->
                        def (fragStart, fragEnd) = obj
                        String dcStatus
                        if (idx == 0) {
                            dcStatus = "C_TERMINAL_DISC"
                        } else if (idx == regions.size() - 1) {
                            dcStatus = "N_TERMINAL_DISC"
                        } else {
                            dcStatus = "NC_TERMINAL_DISC"
                        }
                        fragments.add(new LocationFragment(fragStart, fragEnd, dcStatus))
                    }

                } else {
                    def (fragStart, fragEnd) = regions[0]
                    fragments.add(new LocationFragment(fragStart, fragEnd, "CONTINUOUS"))
                }

                Location location = new Location(start, end, hmmLength, evalue, fragments)
                Match match = matches[seqId][modelId]
                if (match == null) {
                    match = new Match(modelId)
                    match.addLocation(location)
                    match.signature = new Signature(superfamilyAccession)
                    matches[seqId][modelId] = match
                } else {
                    match.addLocation(location)
                }
            }
        }
    }

    def json = JsonOutput.toJson(matches)
    def outputFilePath = task.workDir.resolve("superfamily.json")
    new File(outputFilePath.toString()).write(json)
}
