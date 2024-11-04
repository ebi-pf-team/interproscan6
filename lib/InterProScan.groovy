// Class and methods for validating the user inputs

import java.nio.file.*

class InterProScan {
    static final def PARAMS = [
        [
            name: "input",
            required: true,
            metavar: "<FASTA>",
            description: "path to FASTA file of sequences to be analysed."
        ],
        [
            name: "datadir",
            required: true,
            metavar: "<DATA-DIR>",
            description: "path to data directory."
        ],
        [
            name: "applications",
            metavar: "<APPLICATIONS>",
            description: "comma-separated applications to scan the sequences with. Default: all."
        ],
        [
            name: "formats",
            metavar: "<FORMATS>",
            description: "comma-separated output formats. Available: JSON,TSV,XML. Default: JSON,TSV,XML."
        ],
        [
            name: "outdir",
            metavar: "<OUTDIR>",
            description: "output directory where results will be saved. Default: current working directory."
        ],
        [
            name: "disable-precalc",
            description: "treat input sequences as DNA/RNA."
        ],
        [
            name: "precalc-url",
            metavar: "<URL>",
            description: "URL of the match lookup service. Default: https://www.ebi.ac.uk/interpro/match-lookup."
        ],
        [
            name: "nucleic",
            description: "do no retrieve pre-calculated matches from the match lookup service."
        ],
        [
            name: "goterms",
            description: "include Gene Ontology (GO) mapping in output files."
        ],
        [
            name: "pathways",
            description: "include pathway mapping in output files."
        ],
        [
            name: "help",
            description: "print the help message and exit."
        ],
        // No description -> not displayed in the help message
        [
            name: "batch-size",
            description: null
        ],
        [
            name: "signalp-mode",
            description: null
        ],
        [
            name: "signalp-gpu",
            description: null
        ],
        [
            name: "apps-config",
            description: null
        ],
        [
            name: "x-refs-config",
            description: null
        ],
    ]

    static void validateParams(params, log) {
        def allowedParams = this.PARAMS.collect { it.name.toLowerCase() }

        // Check that all params are recognized
        for (e in params) {
            def paramName = e.key
            def paramValue = e.value

            if (paramName.contains("-")) {
                /*
                    From https://www.nextflow.io/docs/latest/cli.html#pipeline-parameters
                    When the parameter name is formatted using `camelCase`,
                    a second parameter is created with the same
                    value using kebab-case, and vice versa.

                    However, we don't want to evalue the `kebab-case` params.
                    And they will eventually be ignored by NF directly,
                    see https://github.com/nextflow-io/nextflow/pull/4702.
                */
                continue
            }

            // Convert to kebab-case
            def kebabParamName = this.camelToKebab(paramName)
            if (!allowedParams.contains(kebabParamName.toLowerCase())) {
                log.warn "Unrecognised option: --${kebabParamName}. Try '--help' for more information."
            }
        }

        // Check that required params (--input, --datadir) are provided
        this.PARAMS.findAll{ it.required }.each { param ->
            def paramName = kebabToCamel(param.name)
            def paramValue = params[paramName]

            if (paramValue == null) {
                log.error "'--${param.name} ${param.metavar}' is mandatory."
                System.exit(1)
            }
        }
    }

    static String resolveFile(String filePath) {
        Path path = Paths.get(filePath)
        return Files.isRegularFile(path) ? path.toRealPath() : null
    }

    static resolveDirectory(String dirPath, boolean mustExist = false, boolean mustBeWritable = false) {
        Path path = Paths.get(dirPath)

        if (Files.exists(path)) {
            if (!Files.isDirectory(path)) {
                return [null, "Not a directory: ${dirPath}."]
            } else if (mustBeWritable && !Files.isWritable(path)) {
                return [null, "Directory not writable: ${dirPath}."]
            }
            return [path.toRealPath(), null]
        } else if (mustExist) {
            return [null, "Not a directory: ${dirPath}."]
        } else {
            try {
                Files.createDirectories(path)
                return [path.toRealPath(), null]
            } catch (IOException) {
                return [null, "Cannot create directory: ${dirPath}."]
            }
        }
    }

    static validateApplications(String applications, Map appsConfig) {
        if (!applications) {
            // Run all applications
            def appsToRun = appsConfig.findAll{ it ->
                !(it.value.disabled)
            }.keySet().toList()
            return [appsToRun, null]
        }

        // Make a collection of recognized application names
        def allApps = [:]
        appsConfig.each { label, appl ->
            allApps[label] = label
            def stdName = appl.name.toLowerCase().replaceAll("[- ]", "")
            allApps[stdName] = label
            (appl.aliases ?: []).each { alias ->
                def stdAlias = alias.toLowerCase().replaceAll("[- ]", "")
                allApps[stdAlias] = label
            }
        }

        def appsToRun = []
        def appsParam = applications.replaceAll("[- ]", "").split(",").collect { it.trim() }.toSet()
        for (appName in appsParam) {
            def key = appName.toLowerCase()
            if (allApps.containsKey(key)) {
                appsToRun.add(allApps[key])
            } else {
                def error = "Unrecognised application: '${appName}'. Try '--help' to list available applications."
                return [null, error]
            }
        }

        return [appsToRun.toSet().toList(), null]
    }

    static List<String> validateSignalpMode(String signalpMode) {
        if (signalpMode.toLowerCase() !in ['fast', 'slow', 'slow-sequential']) {
            def error = "Unrecognised SignalP mode: '${signalpMode}'. Accepted modes: 'fast', 'slow', 'slow-sequential'"
            return [null, error]
        }
        else {
            return [signalpMode.toLowerCase(), null]
        }
    }

    static String kebabToCamel(String kebabName) {
        return kebabName.split('-').toList().indexed().collect { index, word ->
            index == 0 ? word : word.capitalize()
        }.join('')
    }

    static String camelToKebab(String camelName) {
        return camelName.replaceAll(/([a-z])([A-Z])/, '$1-$2').toLowerCase()
    }

    static void printHelp(appsConfig) {
        def result = new StringBuilder()
        result << "Usage: nextflow run ebi-pf-team/interproscan6 -profile <PROFILE> --input <FASTA> --datadir <DATADIR> \n\n"
        result << "Mandatory parameters:\n"
        result << "  -profile <PROFILE>: use this parameter to choose a configuration profile.\n"

        this.PARAMS.findAll{ it.required }.each { param ->
            result << this.formatOption(param) << "\n"
        }

        result << "\nOptional parameters:\n"
        this.PARAMS.findAll{ !it.required && it.description }.each { param ->
            result << this.formatOption(param) << "\n"
        }

        result << "\nAvailable applications:\n"
        appsConfig.each { label, appl ->
            result << "  ${appl.name.replace(' ', '-')}\n"
        }

        print result.toString()
    }

    static String formatOption(option) {
        def text = "  --${option.name}"
        if (option.metavar) {
            text += " ${option.metavar}"
        }

        return text.padRight(40) + ": ${option.description}"
    }
}
