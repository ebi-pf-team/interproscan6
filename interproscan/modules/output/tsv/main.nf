import groovy.json.JsonSlurper
import java.time.LocalDate

process WRITE_TSV_OUTPUT {
    label 'write_output'

    input:
    val matches
    val outputPath

    exec:
    def outputFilePath = "${outputPath}.ips6.tsv"
    def tsvFile = new File(outputFilePath.toString())
    tsvFile.text = "" // clear the file if it already exists

    // Each line contains: seqId md5 seqLength memberDb modelAcc entryDesc start end evalue status date entryAcc entryDesc xrefs
    def currentDate = LocalDate.now().toString()
    def jsonSlurper = new JsonSlurper()
    jsonSlurper.parse(matches).each { seqData ->
        seqData["matches"].each { modelAccession, matchData ->
            Match match = Match.fromMap(matchData)
            String memberDb = match.signature.signatureLibraryRelease.library
            String goterms = match.signature.entry?.goXRefs ? (match.signature.entry.goXRefs.isEmpty() ? '-' : match.signature.entry.goXRefs.join('|')) : '-'
            String pathways = match.signature.entry?.pathwayXRefs ? (match.signature.entry.pathwayXRefs.isEmpty() ? '-' : match.signature.entry.pathwayXRefs.join('|')) : '-'
            String entryAcc = match.signature.entry?.accession ?: '-'
            String entryDesc = match.signature.entry?.description ?: '-'
            char status = 'T'

            seqData["xref"].each { xrefData ->
                match.locations.each { location ->
                    int start = location.start
                    int end = location.end
                    double evalue = location.evalue
                    switch (memberDb) {
                        case ["prints", "signalp", "signalp_euk"]:
                            evalue = location.pvalue
                        case ["cdd", "hamap", "prosite_profiles"]:
                            evalue = location.score
                        case ["coils", "mobidblite", "phobius", "prosite_patterns"]:
                            evalue = "-"
                        default:
                            evalue = location.evalue
                    }

                    tsvFile.append([
                        xrefData["id"], seqData["md5"], seqData["sequenceLength"],
                        memberDb, entryDesc,
                        start, end, evalue, status,
                        currentDate,
                        entryAcc, entryDesc, goterms, "${pathways}\n"
                    ].join('\t'))
                }
            }
        }
    }
}
