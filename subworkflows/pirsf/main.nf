include { SEARCH_PIRSF; PARSE_PIRSF } from  "../../modules/pirsf"

workflow PIRSF {
    take:
    fasta    // [meta, fasta]
    pirsf    // [hmm, dat]

    main:
    ch_hmm = pirsf.map { hmm, _dat -> hmm }
    ch_dat = pirsf.map { _hmm, dat -> dat }

    ch_parsed = ch_dat.map { dat ->
        def parsed = parseDatFile(dat)
        tuple(parsed[0], parsed[1])
    }

    SEARCH_PIRSF(
        fasta,
        ch_hmm
    )

    PARSE_PIRSF(
        SEARCH_PIRSF.out,
        ch_parsed
    )

    emit:
    PARSE_PIRSF.out  // [ meta, json ]
}

def parseDatFile(datFile) {
    def models = [:]    // PIRSF -> list of 5 doubles (meanL, stdL, minS, meanS, stdS)
    def subfamilies = [:]   // child PIRSF -> parent PIRSF
    new File(datFile.toString()).withReader { reader ->
        def accession = null
        reader.eachLine { line ->
            if (line.startsWith('>')) {
                def parts = line.split(/\s+/)
                accession = parts[0].replace('>', '')

                def match = (line =~ /^>PIRSF\d+\schild:\s(.+)$/)
                if (match) {
                    match[0][1].trim().split(/\s+/).each { child ->
                        subfamilies[child] = accession
                    }
                }

            } else if (line ==~ /\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*\s+\d+\.?\d*/) {
                def values = line.split(/\s+/)*.toDouble()
                models[accession] = values
            }
        }
    }

    return [models, subfamilies]
}