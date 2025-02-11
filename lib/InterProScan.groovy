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
            name: "datadir",  // only required when using members with datafiles
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
        [
            name: "max-workers",
            description: "define maximum number of workers available for the InterProScan when running locally."
        ],
        [
            name: "lookup-host",
            description: "define host to lookup service, just change in case of using a local MLS database."
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
        [
            name: "lookup-service",
            description: null
        ],
        [
            name: "api-chunk-size",
            description: null
        ],
        [
            name: "max-retries",
            description: null
        ],
    ]

    static final def VALID_FORMATS = ["JSON", "TSV", "XML"]

    static final def LICENSED_SOFTWARE = ["phobius", "signalp_euk", "signalp_prok", "deeptmhmm"]

    static final def DATA_TYPE = [
            "FILE": ["cla", "clan", "dat", "disc_regs", "evaluator", "hierarchy", "hmm", "hmmbin",
                     "model", "model2sfs", "pdbj95d", "rules", "seed", "selfhits", "site_annotations",
                     "skip_flagged_profiles"],
            "DIR": ["dir", "msf", "paint", "rpsblast_db", "rpsproc_db"]
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
                log.warn "Unrecognised option: '--${paramName}'. Try '--help' for more information."
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
        if (!dirPath && mustExist) { // triggered when data dir is needed but --datadir not used
            return [null, "'--datadir <DATA-DIR>' is required for the selected applications."]
        }
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
            // Run all applications, except licensed packages with an unpopulated dir field
            def appsToRun = appsConfig.findAll{ it ->
                if (this.LICENSED_SOFTWARE.contains(it.key)) {
                    return it.value?.dir
                }
                return true
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

    static validateAppData(List<String> appsToRun, Path datadir, Map appsConfig) {
        def errorMsg = appsToRun.collectMany { appName ->
            appsConfig[appName].collect { key, value ->
                if (this.DATA_TYPE["FILE"].contains(key)) {
                    if (!resolveFile(datadir.resolve(value).toString())) {
                        return "${appName}: file: '${key}': ${value ?: 'null'}"
                    }
                } else if (this.DATA_TYPE["DIR"].contains(key)) {
                    if (!value) {
                        return "${appName}: dir: '${key}': 'null'"
                    }
                    Path dirPath = this.LICENSED_SOFTWARE.contains(appName) ? Paths.get(value) : datadir.resolve(value)
                    if (!Files.exists(dirPath) || !Files.isDirectory(dirPath)) {
                        return "${appName}: dir: '${key}': ${value ?: 'null'}"
                    }
                }
                return null
            }.findAll { it }
        }.join('\n')
        return errorMsg ? "Could not find the following data files\n${errorMsg}" : null
    }

    static validateXrefFiles(Path datadir, Map xRefsConfig, boolean goterms, boolean pathways) {
        def errorMsg = []
        def addError = { type, suffix ->
            String path = datadir.resolve("${xRefsConfig[type]}${suffix}")
            if (!resolveFile(path)) {
                errorMsg << "${type}${suffix}: ${path}"
            }
        }
        addError('entries', '')  // we hard code the file ext in xrefsconfig so no suffix needed here
        if (goterms) {
            addError('goterms', '.ipr.json')
            addError('goterms', '.json')
        }
        if (pathways) {
            addError('pathways', '.ipr.json')
            addError('pathways', '.json')
        }
        return errorMsg ? "Could not find the following XREF data files\n${errorMsg.join('\n')}" : null
    }

    static Set<String> validateFormats(String userFormats) {
        Set<String> formats = userFormats.toUpperCase().split(',') as Set
        def invalidFormats = formats - VALID_FORMATS
        return invalidFormats ? [null, "Invalid output file format provided:\n${invalidFormats.join('\n')}"] : [formats, null]
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
