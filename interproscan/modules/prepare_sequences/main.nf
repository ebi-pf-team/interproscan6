process VALIDATE_FASTA {
    // check the formating of the intput FASTA, i.e. look for illegal characters
    label         'local', 'ips6_container'
    errorStrategy 'terminate'

    input:
    path fasta
    val isNucleic
    val appsToRun
    val appsConfig

    output:
    path fasta
    path "invalid_char_summary"

    script:
    def commands = ""
    def checkFile = "invalid_char_check"

    // General illegal character check
    def alphabet = isNucleic ? FastaFile.NUCLEIC_ALPHABET : FastaFile.PROTEIN_ALPHABET
    def alphabetName = isNucleic ? "nucleic acid" : "protein"
    commands += "echo 'Check: General: Tolerated characters: $alphabet' >> $checkFile\n"
    commands += "grep -vi '^>' $fasta | grep -Eni '[^$alphabet]' >> $checkFile || true\n"

    // Member database illegal specific character checks
    appsToRun.each { app ->
        def forbiddenChars = appsConfig[app]?.invalid_chars ?: []
        forbiddenChars.each { forbiddenChar ->
            forbiddenChar = forbiddenChar.toString()
            commands += "echo 'Check: $app $forbiddenChar' >> $checkFile\n"
            commands += "grep -Eni '^[^>].*[$forbiddenChar${forbiddenChar.toLowerCase()}]' $fasta >> $checkFile || true\n"
        }
    }

    """
    # Run inside the script block to ensure it is run inside the ips6 container which contains grep
    ${commands}

    touch invalid_char_summary

    # Check if the file has only the "Check" lines (i.e., no grep results, meaning no invalid characters)
    line_count=\$(wc -l < $checkFile)
    check_line_count=\$(grep -c "^Check:" $checkFile)

    if [ "\$line_count" -eq "\$check_line_count" ]; then
        echo "No invalid chars"
        # Only contains the initial Check lines, no errors found
        # Do nothing
        :
    else
        echo "Parsing $checkFile"
        # Build a summary of the invalid characters that were found
        current_app=""
        current_char=""

        while IFS= read -r line; do
            if [[ \$line == Check:* ]]; then
                # Extract app info - for general check or specific app check
                if [[ \$line == *"Tolerated characters:"* ]]; then
                    current_app="general"
                    current_char=""
                else
                    # For app-specific checks like "Check: antifam -"
                    current_app=\$(echo "\$line" | awk '{print \$2}')
                    current_char=\$(echo "\$line" | awk '{print \$3}')
                fi

                # Reset line number collection for new check
                line_nums=""
            elif [[ \$line =~ ^[0-9]+: ]]; then
                # Extract line number from grep output
                line_num=\$(echo "\$line" | grep -o '^[0-9]*')

                # Append to line numbers for this check
                if [ -z "\$line_nums" ]; then
                    line_nums="\$line_num"
                else
                    line_nums="\$line_nums, \$line_num"
                fi

                # Write to summary when we have line numbers
                if [ -n "\$current_app" ] && [ -n "\$line_nums" ]; then
                    if [ -n "\$current_char" ]; then
                        echo "\$current_app forbidden character '\$current_char' found on lines: \$line_nums" >> invalid_char_summary
                    else
                        echo "\$current_app forbidden character(s) found on lines: \$line_nums" >> invalid_char_summary
                    fi
                fi
            fi
        done < $checkFile
    fi
    """
}

process LOAD_SEQUENCES {
    // Populate a local sqlite3 database with sequences from the pipeline's input FASTA file.
    label         'local'
    errorStrategy 'terminate'

    input:
    val fasta
    val nucleic

    output:
    path "sequences.db"

    exec:
    def outputFilePath = task.workDir.resolve("sequences.db")
    SeqDB db = new SeqDB(outputFilePath.toString())
    db.loadFastaFile(fasta.toString(), nucleic, false)
    db.close()
}

process LOAD_ORFS {
    // add protein seqs translated from ORFS in the nt seqs to the database
    label         'local'
    errorStrategy 'terminate'

    input:
    val translatedFastas  // could be one or multiple paths
    val dbPath

    output:
    val dbPath // ensure BUILD_BATCHES runs after LOAD_ORFS

    exec:
    SeqDB db = new SeqDB(dbPath.toString())
    translatedFastas.each {
        db.loadFastaFile(it.toString(), false, true)
    }
    db.close()
}

process SPLIT_FASTA {
    // Build the FASTA file batches of unique protein sequences for the sequence analysis
    label         'local'
    errorStrategy 'terminate'

    input:
    val dbPath
    val batchSize
    val nucleic

    output:
    path "*.fasta"

    exec:
    String prefix = task.workDir.resolve("input").toString()
    SeqDB db = new SeqDB(dbPath.toString())
    db.splitFasta(prefix, batchSize, nucleic)
    db.close()
}
