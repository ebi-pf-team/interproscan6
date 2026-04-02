package uk.ac.ebi.interpro

import groovy.json.JsonSlurper
import java.security.MessageDigest
import java.nio.file.*
import uk.ac.ebi.interpro.HTTPRequest

class InterProScan {
    static final def PARAMS = [
        [
            name: "applications",
            metavar: "<LIST>",
            description: "Comma-separated list of applications to run. Default: all available applications.",
        ],
        [
            name: "cpus",
            metavar: "<N>",
            description: "Number of CPU threads allocated per task for applications that support multithreading. Default: 1."
        ],
        [
            name: "datadir",  // only required when using members with datafiles
            metavar: "<DATADIR>",
            description: "Path to the data directory. Required when running applications that depend on external data files (e.g. Pfam, CDD). Not required when using only self-contained applications (e.g. COILS, MobiDB-lite).",
        ],
        [
            name: "formats",
            metavar: "<LIST>",
            description: "Comma-separated output formats. Available: gff3, json, jsonl, tsv, xml. Default: all formats."
        ],
        [
            name: "globus",
            description: "Use the Globus mirror of the EMBL-EBI FTP site to download InterProScan data files.",
        ],
        [
            name: "goterms",
            description: "Include Gene Ontology (GO) annotations in output files."
        ],
        [
            name: "help",
            description: "Output the usage message and exit."
        ],
        [
            name: "input",
            required: true,
            metavar: "<FASTA>",
            description: "Path to FASTA file of sequences to be analysed."
        ],
        [
            name : "interpro",
            metavar: "<VERSION>",
            description: "InterPro release version to use. Default: latest available."
        ],
        [
            name: "matches-api-url",
            metavar: "<URL>",
            description: "URL of an alternative InterPro Matches API instance to retrieve precomputed results.",
        ],
        [
            name: "max-workers",
            metavar: "<N>",
            description: "Maximum number of parallel tasks executed when running locally."
        ],
        [
            name: "no-matches-api",
            description: "Disable use of the Matches API. All analyses are executed locally."
        ],
        [
            name: "nucleic",
            description: "Interpret input sequences as nucleotides. Sequences are translated in all six reading frames and resulting ORFs are analysed."
        ],
        [
            name: "outdir",
            metavar: "<DIR>",
            description: "Directory where output files are written. Default: current working directory."
        ],
        [
            name: "outprefix",
            metavar: "<PREFIX>",
            description: "Base name for output files (no directory or extension). Default: input filename."
        ],
        [
            name: "pathways",
            description: "Include pathway annotations in output files."
        ],
        [
            name: "run-ml",
            description: "Enable machine learning (ML) based analyses (e.g. DeepTMHMM, SignalP, TMbed). Disabled by default due to high resource usage."
        ],
        [
            name: "skip-applications",
            metavar: "<LIST>",
            description: "Comma-separated list of applications to exclude from the analysis. Default: none.",
        ],
        [
            name: "skip-interpro-version-check",
            description: "Disable the version compatibility check between InterProScan and InterPro data."
        ],
        [
            name: "use-gpu",
            description: "Enable GPU acceleration for supported ML-based applications."
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
            name: "sub-batch-size",
            description: null
        ],
        [
            name: "skip-repr-locations",
            description: null
            // Used in production. Skips identifying representative locations
        ],
        [
            name: "apps-config",
            description: null
        ],
        [
            name: "matches-api-chunk-size",
            description: null
        ],
        [
            name: "matches-api-max-retries",
            description: null
        ],
    ]

    static final def VALID_FORMATS = ["JSON", "JSONL", "TSV", "XML", "GFF3"]

    static final def LICENSED_SOFTWARE = ["interpro_n", "phobius", "signalp_euk", "signalp_prok", "deeptmhmm"]

    static final def DATA_TYPE = [
            "FILE": ["cla", "clan", "dat", "disc_regs", "evaluator", "hierarchy", "hmm", "hmmbin",
                     "model", "model2sfs", "pdbj95d", "rules", "seed", "selfhits", "site_annotations",
                     "skip_flagged_profiles"],
            "DIR": ["msf", "paint", "rpsblast_db", "rpsproc_db"]
    ]

    static final String FTP_URL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"
    static final String GLOBUS_URL = "https://g-a8b222.dd271.03c0.data.globus.org/pub/software/unix/iprscan/6"

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
                if (paramObj?.metavar != null && paramValue != null && (paramValue instanceof Boolean)) {
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

    static String getMD5Hash(Path file) {
        // Get the MD5 Hash of a local file. Used to check if the correct or complete file has been downloaded
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

    static List<String> getAppsWithData(List<String> applications, Map appsConfig) {
        return applications.findAll { String appName ->
            appsConfig.get(appName)?.has_data
        }
    }

    static findLocalHighestVersionDir(String dirPath) {
        Path path = Paths.get(dirPath)
        if (!Files.exists(path)) {  // check if exists, otherwise it will raise a generic "No such file or dir" err
            return null
        }
        def dirs = Files.list(path)
            .findAll { Files.isDirectory(it) && it.fileName.toString() ==~ /^\d+\.\d+$/ }
            .sort { a, b ->
                def (aMajor, aMinor) = a.fileName.toString().tokenize('.').collect { it.toInteger() }
                def (bMajor, bMinor) = b.fileName.toString().tokenize('.').collect { it.toInteger() }
                return aMajor <=> bMajor ?: aMinor <=> bMinor
            }
        return dirs ? dirs.last().fileName.toString() : null
    }

    static List<String>validateApplications(String applications, String skipApplications, Map appsConfig, Boolean runML) {
        // Return: [appsToRun | null, errorMessage | null, warningMessage | null]
        if (applications && runML) {
            return [null, "--applications and --run-ml are mutually exclusive", null]
        }
        if (applications && skipApplications) {
            return [null, "--applications and --skip-applications are mutually exclusive", null]
        }

        // Normalize "catalogue" of applications
        Map<String, String> aliases = [:]
        appsConfig.each { label, cfg ->
            def std = cfg.name.toLowerCase().replaceAll(/[-_ ]/, "")
            aliases[std] = label
            (cfg.aliases ?: []).each { alias ->
                def stdAlias = alias.toLowerCase().replaceAll(/[-_ ]/, "")
                aliases[stdAlias] = label
            }
        }

        Set<String> mlAppls = appsConfig.findAll { k, v ->
            v.containsKey("enabled")
        }.keySet() as Set

        Closure<Boolean> isAvailable = { label ->
            def cfg = appsConfig[label]
            assert cfg != null
            if (this.LICENSED_SOFTWARE.contains(label)) {
                return cfg.dir != null && !cfg.dir.isEmpty()
            }
            return true
        }

        // Initialize selection
        Set<String> base = [] as Set
        if (!applications) {
            // default: include all non-ML, and add ML if runML is set
            appsConfig.keySet().each { label ->
                boolean isML = mlAppls.contains(label)
                if (isAvailable(label) && (!isML || (isML && runML))) {
                    base << label
                }
            }
        }

        // Explicit --applications
        if (applications) {
            Set<String> req = applications.split(",")
                                          .collect { it.trim().toLowerCase().replaceAll(/[-_ ]/, "") }
                                          .findAll { it }
                                          .toSet()
            Set<String> unrec = []
            Set<String> unavail = []
            req.each { raw ->
                def label = aliases[raw]
                if (!label) {
                    unrec << raw
                    return
                }
                if (!isAvailable(label)) {
                    unavail << raw
                    return
                }
                base << label
            }

            if (!unrec.isEmpty() || !unavail.isEmpty()) {
                def msg = []
                if (!unrec.isEmpty())
                    msg << "Unrecognised applications: ${unrec.join(', ')}"
                if (!unavail.isEmpty())
                    msg << "Unavailable applications: ${unavail.join(', ')}"
                return [null, msg.join('. '), null]
            }
        }
        
        // Explicit --skip--applications
        if (skipApplications) {
            Set<String> req = skipApplications.split(",")
                                              .collect { it.trim().toLowerCase().replaceAll(/[-_ ]/, "") }
                                              .findAll { it }
                                              .toSet()
            Set<String> unrec = []
            req.each { raw ->
                def label = aliases[raw]
                if (!label) {
                    unrec << raw
                    return
                }
                base.remove(label)
            }

            if (!unrec.isEmpty()) {
                def msg = "Unrecognised applications: ${unrec.join(', ')}"
                return [null, msg, null]
            }
        }

        // Warnings
        if (runML && !applications) {
            // ML requested; warn if any ML is unavailable
            def unavailableML = mlAppls.findAll { label ->
                !isAvailable(label)
            }
            if (!unavailableML.isEmpty()) {
                def warn = "Unavailable ML analyses will be skipped: ${unavailableML.join(', ')}"
                return [base.toList(), null, warn]
            }
        }

        return [base.toList(), null, null]
    }

    static enableGpuAcceleration(Boolean useGpu, List<String> applications, Map appsConfig) {
        // return apps_config, warn
        if (!useGpu) {
            return [appsConfig, null]
        }

        // no gpu compatible apps requested
        def gpuApps = applications.findAll { appName ->
            appsConfig[appName]?.containsKey('use_gpu')
        }
        if (!gpuApps) {
            return [appsConfig, "'--use-gpu' was specified but no GPU-compatible applications were requested, so GPU acceleration will not be used."]
        }
        
        // enable gpu for requested apps
        def updatedAppsConfig = appsConfig.collectEntries { appName, thisAppConfig ->
            if (gpuApps.contains(appName)) {
                thisAppConfig.use_gpu = true
            }
            return [(appName): thisAppConfig]
        }
        return [updatedAppsConfig, null]
    }

    static String validateInterProVersion(versionParam) {
        String version = null
        if (versionParam instanceof Number) {
            double num = ((Number) versionParam).doubleValue();
            if (num == Math.floor(num)) {
                version = String.format("%.1f", num);
            } else {
                version = Double.toString(num);
            }
        } else if (versionParam instanceof String && versionParam == "latest") {
            version = versionParam
        }
        return version
    }

    static List<String> fetchCompatibleVersions(String majorMinorVersion, boolean useGlobus = false) {
        String baseUrl = useGlobus ? GLOBUS_URL : FTP_URL
        String url = "${baseUrl}/${majorMinorVersion}/versions.json"
        Map versions = HTTPRequest.fetch(url, null, 2, false)
        return versions?.interpro?.collect { it?.toString() } ?: null
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
        def optionsMap = this.PARAMS
            .findAll { it.description }
            .collectEntries { opt ->
                [
                    (opt.name): [
                        description: opt.description,
                        metavar    : opt.metavar
                    ]
                ]
            }

        def result = new StringBuilder()
        result << "Usage: nextflow run ebi-pf-team/interproscan6 -profile <PROFILE> --input <FASTA> [--datadir <DATADIR>] [options]\n\n"
        result << "Mandatory parameters:\n"
        result << this.formatOption("profile",
                                    [
                                        metavar    : "<PROFILE>", 
                                        description: "Execution profile defining the runtime environment (e.g. slurm, lsf, docker)."
                                    ],
                                    " -")
        result << this.formatOption("input", optionsMap.input)

        result << "\nConditional parameters:\n"
        ["datadir"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nExecution and resource control:\n"
        ["cpus", "max-workers", "use-gpu"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nInput interpretation:\n"
        ["nucleic"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nAnalysis selection:\n"
        ["applications", "skip-applications", "run-ml"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nOutput control:\n"
        ["outdir", "outprefix", "formats", "goterms", "pathways"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nExternal services and data:\n"
        ["interpro", "matches-api-url", "no-matches-api", "globus"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }

        result << "\nOther options:\n"
        ["help"].each { paramName ->
            assert optionsMap.containsKey(paramName)
            result << this.formatOption(paramName, optionsMap[paramName])
        }
        
        result << "\nAvailable applications:\n"
        appsConfig.each { label, appl ->
            if (appl.name != "InterPro-N") {
                result << "  ${appl.name.replace(' ', '-')}\n"
            }
        }

        print result.toString()
    }

    static String formatOption(String name, Map<String, String> props, String prefix = "--") {
        def text = "  ${prefix}${name}"
        if (props.metavar) {
            text += " ${props.metavar}"
        }

        assert props.description != null
        return text.padRight(35) + "${props.description}\n"
    }

    static List parseAppsConfig(Boolean useGpu, List<String> apps, Path appsConfigFile) {
        ConfigSlurper configSlurper = new ConfigSlurper()
        def config = configSlurper.parse(appsConfigFile.toURI().toURL())
        def warn = null
        def appsConfig = config.params.appsConfig
        
        // set all keys to lowercase to simplify matching with user input
        appsConfig = appsConfig.collectEntries { appName, cfg ->
            [(appName.toLowerCase()): cfg]
        }
        
        (appsConfig, warn) = this.enableGpuAcceleration(useGpu, apps, appsConfig)
        return [appsConfig, warn]
    }
}
