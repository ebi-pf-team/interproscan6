// Class and methods for validating the user inputs

import java.security.MessageDigest
import java.nio.file.*
import InterPro

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
            description: "path to data directory.",
            canBeNull: true
        ],
        [
            name: "applications",
            metavar: "<APPLICATIONS>",
            description: "comma-separated applications to scan the sequences with. Default: all.",
            canBeNull: true
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
            name: "offline",
            description: "run InterProScan in offline mode, disabling queries to the InterPro Matches API. Pre-calculated matches for known sequences will not be retrieved, and analyses will be run locally."
        ],
        [
            name : "interpro",
            description: "the InterPro release to be used. Defaults to 'latest'."
        ],
        [
            name: "matches-api-url",
            metavar: "<URL>",
            description: "override the default InterPro Matches API, hosted at EMBL-EBI. Use this option to specify the URL of an alternative Matches API instance.",
            canBeNull: true
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
            name: "download",
            description: ("Download data for the selected applications. " +
                    "If `<DATA-DIR>/xrefs` exists, InterProScan will use the existing InterPro data. " +
                    "Otherwise, the latest compatible InterPro release will be downloaded."
            )
        ],
        [
            name: "help",
            description: "print the help message and exit."
        ],
        [
            name: "max-workers",
            description: "define maximum number of workers available for the InterProScan when running locally."
        ],
        /*
        If an option's description is set to null, it will be hidden from the help message
        and no "Unrecognised option" warning will be produced.
        Use this for params defined in config files that should not be available on the command line
        */
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
            name: "ftp-server",
            description: null
        ],
        [
            name: "lookup-service",
            description: null
        ],
    ]

    static final def VALID_FORMATS = ["JSON", "TSV", "XML"]

    static final def LICENSED_SOFTWARE = ["phobius", "signalp_euk", "signalp_prok", "deeptmhmm"]

    static final def DATA_TYPE = [
            "FILE": ["cla", "clan", "dat", "disc_regs", "evaluator", "hierarchy", "hmm", "hmmbin",
                     "model", "model2sfs", "pdbj95d", "rules", "seed", "selfhits", "site_annotations",
                     "skip_flagged_profiles"],
            "DIR": ["msf", "paint", "rpsblast_db", "rpsproc_db"]
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
            if (allowedParams.contains(kebabParamName.toLowerCase())) {
                def paramObj = this.PARAMS.find { it.name.toLowerCase() == kebabParamName.toLowerCase() }
                assert paramObj != null
                if (paramObj?.metavar != null && !paramObj?.canBeNull && !(paramValue instanceof String)) {
                    log.error "'--${paramObj.name} ${paramObj.metavar}' is mandatory and cannot be empty."
                    System.exit(1)
                }
            } else {
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

    static String getMD5Hash(String filePath) {
        // Get the MD5 Hash of a local file. Used to check if the correct or complete file has been downloaded
        def file = new File(filePath)
        MessageDigest md = MessageDigest.getInstance("MD5")
        file.withInputStream { is ->
            byte[] buffer = new byte[8192]
            int bytesRead
            while ((bytesRead = is.read(buffer)) != -1) {
                md.update(buffer, 0, bytesRead)
            }
        }
        return md.digest().collect { String.format("%02x", it)}.join()
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

    static validateAppData(List<String> appsToRun, Path datadir, Map appsConfig, Map databaseMap, Boolean returnList=false) {
        def missingApps = [] as Set // only returned if returnList is true
        def errorMsg = appsToRun.collectMany { appName ->
            def fmtAppName = InterPro.formatMemberDbName(appName)

            def appVersion = databaseMap[fmtAppName].toString()

            if (!appVersion) {
                return "Could not retrieve release for $appName"
            }
            if (!appsConfig[appName]['dir']) {  // member does not have a datadir, e.g. coils and mobidblite
                return null
            }
            def appDir = datadir.resolve("${appsConfig[appName]['dir']}/$appVersion".toString())
            if (!Files.exists(appDir) || !Files.isDirectory(appDir)) {
                missingApps.add(appName)
                return ["${appName}: dir ${appDir.toString()} is missing"]
            }
            appsConfig[appName].collect { key, value ->
                if (this.DATA_TYPE["FILE"].contains(key)) {
                    if (!resolveFile(appDir.resolve(value).toString())) {
                        missingApps.add(appName)
                        return "${appName}: file: '${key}': ${value ?: 'null'}"
                    }
                } else if (this.DATA_TYPE["DIR"].contains(key)) {
                    if (!value) {
                        missingApps.add(appName)
                        return "${appName}: dir: '${key}': 'null'"
                    }
                    Path dirPath = this.LICENSED_SOFTWARE.contains(appName) ? Paths.get(value) : appDir.resolve(value)
                    if (!Files.exists(dirPath) || !Files.isDirectory(dirPath)) {
                        missingApps.add(appName)
                        return "${appName}: dir: '${key}': ${value ?: 'null'}"
                    }
                }
                return null
            }.findAll { it }
        }.join('\n')
        return returnList ? missingApps as List : (errorMsg ? "Could not find the following data files\n${errorMsg}" : null)
    }

    static validateXrefFiles(Path datadir, String interproRelease, Map xRefsConfig, boolean goterms, boolean pathways) {
        def error = ""
        def addError = { type, suffix ->
            String path = datadir.resolve("interpro/$interproRelease/${xRefsConfig[type]}${suffix}").toString()
            if (!resolveFile(path)) {
                error << "${type}${suffix}: ${path}"
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
        return error ? "Could not find the following XREF data files\n${error.join('\n')}" : null
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
