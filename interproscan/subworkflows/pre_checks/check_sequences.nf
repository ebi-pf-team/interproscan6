def isNucleic(sequence) {
    def LEGAL_NUCLEIC_CHARS = ['A', 'T', 'C', 'G', 'U', '*', '-'] as Set
    sequence.toUpperCase().every { LEGAL_NUCLEIC_CHARS.contains(it) }
}

def checkIllegalChars(sequence, apps, nucleic) {
    def ILLEGAL_CHARS = [
        "antifam": "-",
        "cdd": "",
        "coils": "",
        "funfam": "-_.",
        "gene3d": "-_.",
        "hamap": "-_",
        "mobidb": "",
        "ncbifam": "-",
        "panther": "-",
        "pfam": "-",
        "phobius": "-*._oxuzj",
        "pirsf": "-",
        "pirsr": "-",
        "prints": "-._",
        "prosite_patterns": "",
        "prosite_profiles": "-._",
        "sfld": "-._",
        "smart": "",
        "superfamily": "-",
        "signalp": "",
        "signalp_euk": ""
    ]
    def errors = [:]
    def invalidCharPattern = ~/[^A-Za-z_\-\*\.]*/
    def invalidChars = sequence.findAll(invalidCharPattern).flatten().unique().collect { "'${it}'" }.findAll { it != "''" }
    if (invalidChars) {
        errors['GENERAL'] = new TreeSet(invalidChars)
    }
    apps.each { app ->
        def appChars = sequence.toSet().intersect(ILLEGAL_CHARS[app].toSet()).collect { "'${it}'" }.findAll { it != "''" }
        if (appChars) {
            if (nucleic && app == 'phobius' && "'u'" in appChars) {
                appChars.remove("'u'")
            }
            if (appChars) {
                errors[app.toUpperCase()] = new TreeSet(appChars)
            }
        }
    }
    return errors.isEmpty() ? null : errors
}

workflow CHECK_SEQUENCES {
    take:
    seq_input
    applications_lower
    is_nucleic

    main:
    def errors = [
        'non-nucleic-chars': new TreeSet(),
        'illegal-chars': [:]
    ]

    def currentSeq = null
    seq_input.eachLine { line ->
        if (line.startsWith('>')) {
            if (currentSeq) {
                // Process the previous sequence
                if (is_nucleic) {
                    if (!isNucleic(currentSeq.sequence)) {
                        errors['non-nucleic-chars'] << currentSeq.seqKey
                    }
                }

                def illegalCharsDetected = checkIllegalChars(currentSeq.sequence, applications_lower, is_nucleic)
                if (illegalCharsDetected) {
                    errors['illegal-chars'][currentSeq.seqKey] = illegalCharsDetected
                }
            }

            currentSeq = [seqKey: line.substring(1).trim(), sequence: '']
        } else {
            currentSeq.sequence += line.trim()
        }
    }

    // Process the last sequence
    if (is_nucleic) {
        if (!isNucleic(currentSeq.sequence)) {
            errors['non-nucleic-chars'] << currentSeq.seqKey
        }
    }

    illegalCharsDetected = checkIllegalChars(currentSeq.sequence, applications_lower, is_nucleic)
    if (illegalCharsDetected) {
        errors['illegal-chars'][currentSeq.seqKey] = illegalCharsDetected
    }

    // If any issues were detected then raise and error and terminate
    if (errors['non-nucleic-chars'] || errors['illegal-chars']) {
        def outdir = file(params.outdir)
        if (!outdir.exists()) {
            outdir.mkdirs()
        }

        def outpath = outdir.resolve('illegal_characters.err')

        def errorFile = file(outpath)
        errorFile.text = ""  // Clear the file if it already exists

        if (errors['non-nucleic-chars']) {
            errorFile.append("Non-nucleic acid residues were found in the following sequences:\n")
            errors['non-nucleic-chars'].each { errorFile.append("${it}\n") }
            errorFile.append("\n")
        }
        if (errors['illegal-chars']) {
            errorFile.append("Illegal characters were detected in the following sequences:\n")
            errors['illegal-chars'].each { seqKey, seqErrors ->
                seqErrors.each { app, chars ->
                    def appPlaceholder = app != "GENERAL" ? " for ${app}" : ""
                    errorFile.append("Sequence ${seqKey.split()[0]} contains illegal character(s)${appPlaceholder}: ${chars.join(', ')}\n")
                }
            }
        }

        if (errors['non-nucleic-chars'] && errors['illegal-chars']) {
            exit 22, "Illegal characters were detected and non-nucleic acid residues were found in the input file:\n${seq_input}\nCheck ${outpath} for details."
        } else if (errors['non-nucleic-chars']) {
            exit 22, "Non-nucleic acid residues were found in the input file:\n${seq_input}\nCheck ${outpath} for details."
        } else if (errors['illegal-chars']) {
            exit 22, "Illegal characters were detected in the input file:\n${seq_input}\nCheck ${outpath} for details."
        }
    }
}
