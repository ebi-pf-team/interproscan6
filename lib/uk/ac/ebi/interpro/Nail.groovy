package uk.ac.ebi.interpro

import java.nio.file.Path


class Nail {
    /*
        Parse the alignment output produced by `nail search --ali-out` (the `.ali` file).

        nail reports one hit per query(HMM)/target(sequence) alignment and, unlike HMMER3,
        does not distinguish sequence-level from domain-level statistics. 
        Each hit maps to a single Location and the parent Match inherits that hit's e-value/score/bias
        (the best hit's values when several align the same model to the same sequence).

        An entry looks like:
            query:        An_peroxidase
            target:       sp|Q9ES45|DUOX2_RAT
            query start:  2
            query end:    530
            target start: 36
            target end:   560
            score:        476.2
            comp bias:    0.0
            E-value:      6.6e-145
            cell frac:    0.015

            ==

                  An_peroxidase     2 rtidGscNnlkn...   80
                                      +++dG++Nnlk ...
            sp|Q9ES45|DUOX2_RAT      36 QRYDGWFNNLKY...  113
                                      79************...
            ...
            //

    */
    static parseOutput(Path filePath, String memberDb) {
        SignatureLibraryRelease libraryRelease = new uk.ac.ebi.interpro.SignatureLibraryRelease(memberDb, null)
        def hits = [:].withDefault { [:] }

        // State for the entry currently being read
        def meta = [:]
        StringBuilder queryAlignment = new StringBuilder()
        StringBuilder targetAlignment = new StringBuilder()
        boolean inAlignment = false

        Closure flush = {
            String modelAccession = meta["query"]
            String targetId = meta["target"]
            if (modelAccession == null || targetId == null) {
                return
            }

            Double evalue = meta["E-value"]?.isEmpty() ? null : meta["E-value"] as Double
            Double score = meta["score"]?.isEmpty() ? null : meta["score"] as Double
            Double bias = meta["comp bias"]?.isEmpty() ? null : meta["comp bias"] as Double

            // nail reports a single interval per hit, so it doubles as the envelope.
            // hmmLength and hmmBounds are not reported by nail (left null).
            Location location = new uk.ac.ebi.interpro.Location(
                meta["target start"].toInteger(),  // start (on the target sequence)
                meta["target end"].toInteger(),    // end
                meta["query start"].toInteger(),   // hmmStart (position on the model)
                meta["query end"].toInteger(),     // hmmEnd
                null,                              // hmmLength
                null,                              // hmmBounds
                meta["target start"].toInteger(),  // envelopeStart
                meta["target end"].toInteger(),    // envelopeEnd
                evalue,
                score,
                bias
            )
            if (queryAlignment.length() > 0) {
                location.queryAlignment = queryAlignment.toString()
            }
            if (targetAlignment.length() > 0) {
                location.targetAlignment = targetAlignment.toString()
            }

            Match match = hits[targetId][modelAccession]
            if (match == null) {
                match = new uk.ac.ebi.interpro.Match(
                    modelAccession,
                    evalue,
                    score,
                    bias,
                    new uk.ac.ebi.interpro.Signature(modelAccession, libraryRelease)
                )
                hits[targetId][modelAccession] = match
            } else if (evalue != null && (match.evalue == null || evalue < match.evalue)) {
                // Several domains of the same model on the same sequence: keep the best
                // hit's statistics at the match level.
                match.evalue = evalue
                match.score = score
                match.bias = bias
            }
            match.addLocation(location)
        }

        filePath.withReader { reader ->
            String line
            while ((line = reader.readLine()) != null) {
                String trimmed = line.trim()
                if (trimmed == "//") {
                    // End of the current entry
                    flush()
                    meta.clear()
                    queryAlignment.setLength(0)
                    targetAlignment.setLength(0)
                    inAlignment = false
                } else if (trimmed == "==") {
                    // Separator between the metadata header and the alignment blocks
                    inAlignment = true
                } else if (!inAlignment) {
                    // Metadata header: "key: value" pairs
                    int idx = trimmed.indexOf(":")
                    if (idx > 0) {
                        meta[trimmed.substring(0, idx).trim()] = trimmed.substring(idx + 1).trim()
                    }
                } else if (!trimmed.isEmpty()) {
                    /*
                        Alignment blocks. Each block has four lines:
                            <query name>  <start> <query alignment>  <end>
                                                  <consensus/match annotation>
                            <target name> <start> <target alignment> <end>
                                                  <posterior probability annotation>
                        The query and target lines start with their respective name, which
                        lets us pick them out and skip the two annotation lines.
                    */
                    def fields = trimmed.split(/\s+/)
                    if (fields.length >= 4) {
                        if (fields[0] == meta["query"]) {
                            queryAlignment.append(fields[2])
                        } else if (fields[0] == meta["target"]) {
                            targetAlignment.append(fields[2])
                        }
                    }
                }
            }
            // Flush a trailing entry that is not terminated by "//"
            flush()
        }

        return hits
    }
}
