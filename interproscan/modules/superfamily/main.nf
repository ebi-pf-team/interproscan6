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

    // output:
    // tuple val(meta), path("antifam.json")    

    exec:
    def matches = [:]

    file(superfamily_out.toString()).eachLine { line ->
        line = line.trim()
        if (line) {
            def fields = line.split(/\s+/)
            assert fields.size() == 9
            String seqId = fields[0]
            String modelId = fields[1]
            if (modelId != "-") {
                String regions = fields[2]
                Double evalue = Double.parseDouble(fields[3])
                Integer hmmStart =  Integer.parseInt(fields[4])
                String alignment = fields[5]
                Double familyEvalue = Double.parseDouble(fields[6])
                Integer scopDomainId = Integer.parseInt(fields[7])
                Integer scopFamilyId = Integer.parseInt(fields[8])

                if (seqId == "tr|Q54AJ4" && modelId == "0040488") {
                    def fragments = []
                    regions.split(",").each { region ->
                        def boundaries = region.split("-")
                        assert boundaries.size() == 2
                        int start = boundaries[0].toInteger()
                        int end = boundaries[1].toInteger()
                        println(boundaries)
                        new LocationFragment(start, end, "CONTINUOUS")
                    }
                }

                
            }
        }
    }

    // def json = JsonOutput.toJson(matches)
    // new File(outputFilePath.toString()).write(json)
}