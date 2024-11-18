import groovy.json.JsonSlurper
import java.time.format.DateTimeFormatter
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
    // Each line contains: seqId md5 seqLength memberDb modelAcc sigDesc start end evalue status date entryAcc entryDesc xrefs
    def currentDate = LocalDate.now().format(DateTimeFormatter.ofPattern("dd-MM-yyyy"))
    def jsonSlurper = new JsonSlurper()
    jsonSlurper.parse(matches).each { seqData ->
        seqData["matches"].each { modelAccession, matchData ->
            Match match = Match.fromMap(matchData)
            String memberDb = match.signature.signatureLibraryRelease.library
            String sigDesc = match.signature.description?: '-'
            String goterms = match.signature.entry?.goXRefs ? (match.signature.entry.goXRefs.isEmpty() ? '-' : match.signature.entry.goXRefs.join('|')) : '-'
            String pathways = match.signature.entry?.pathwayXRefs ? (match.signature.entry.pathwayXRefs.isEmpty() ? '-' : match.signature.entry.pathwayXRefs.join('|')) : '-'
            String entryAcc = match.signature.entry?.accession ?: '-'
            String entryDesc = match.signature.entry?.description ?: '-'
            char status = 'T'

            seqData["xref"].each { xrefData ->
                match.locations.each { Location loc ->
                    int start = loc.start
                    int end = loc.end
                    def scoringValue = "-"
                    switch (memberDb) {
                        case ["cdd", "prints"]:
                            scoringValue = match.evalue
                            break
                        case ["signalp", "signalp_euk"]:
                            scoringValue = loc.pvalue
                            break
                        case ["hamap", "prositeprofiles"]:
                            scoringValue = loc.score
                            break
                        case ["coils", "mobidblite", "phobius", "prositepatterns"]:
                            scoringValue = "-"
                            break
                        default:
                            scoringValue = loc.evalue
                            break
                    }

                    tsvFile.append([
                        xrefData["id"], seqData["md5"], seqData["sequence"].length(),
                        Output.convertDbName(memberDb),
                        match.signature.accession, sigDesc,
                        start, end, scoringValue, status,
                        currentDate,
                        entryAcc, entryDesc, goterms, "${pathways}\n"
                    ].join('\t'))
                }
            }
        }
    }
}
