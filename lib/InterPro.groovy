// Class and methods for interacting with InterPro servers and data

import groovy.json.JsonSlurper
import java.io.File
import java.io.IOException
import java.net.URL  // TODO: Depracted. Migrate to URI
import java.nio.file.*
import InterProScan

class InterPro {
    static final def FTP_URL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"

    static Map getMemberDbReleases(Map xRefsConfig, String interproRelease, datadir) {
        /* Load the datadir/interpro/database.json file and set all keys to lowercase to match applications.config
        Don't worry about checking it exists, this was done in InterProScan.validateXrefFiles() */
        JsonSlurper jsonSlurper = new JsonSlurper()
        def databaseJsonPath = datadir.toString() + "/${xRefsConfig.dir}/${interproRelease}/${xRefsConfig.databases}"
        System.out.println( "** $databaseJsonPath -- databaseJsonPath" )
        def databaseJson = new File(databaseJsonPath.toString())
        def memberDbReleases = jsonSlurper.parse(databaseJson)
        memberDbReleases = memberDbReleases.collectEntries { appName, versionNum ->
            [(appName.toLowerCase()): versionNum]
        }
        return memberDbReleases
    }

    static List<String> checkCompatibility(String iprscan, String interpro, def log) {
        // Check that the InterPro and InterProScan versions are compatible
        def noConnWarning = ""
        def latestReleaseWarning = ""
        def error = ""
        def releasesURL = InterPro.sanitizeURL("${FTP_URL}/${iprscan}/releases.json")
        Map compatibleReleases = InterPro.httpRequest(releasesURL, null, 0, true, log)
        if (!compatibleReleases) {
            noConnWarning = "An error occurred while querying the EBI ftp to assess the compatibility of\n"+
                     "InterPro release '${interpro}' with InterProScan version '${iprscan}'\n"+
                     "Proceeding with using InterPro release '${interpro}'.\n"+
                    "Warning, InterProScan may not behave as expected"
        } else if (!compatibleReleases[iprscan].contains(interpro)) {
            error = "InterPro release '$interpro' is not compatiable with "
        } else if (interpro < compatibleReleases[iprscan].max()) {
            latestReleaseWarning = "A later compatible InterPro release '${compatibleReleases[iprscan].max()}' is availble.\n"+
                                   "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release"
        }
        return [noConnWarning, latestReleaseWarning, error]
    }

    static formatMemberDbName(String memberName) {
        def fmtMemberName = memberName.toLowerCase()
        if (fmtMemberName.startsWith("cath")) {
            return fmtMemberName.replace("cath", "cath-")
        } else if (fmtMemberName.startsWith("prosite")) {
            return  fmtMemberName.replace("prosite", "prosite ")
        } else {
            return fmtMemberName
        }
    }

    static getInterproRelease(Path datadir, String interpro) {
        // Checks the InterProDir exists and then retrieve the specified interpro release number
        def _interproRelease = null
        def error = ""

        // Check if the <datadir>/interpro dir exists
        def _interproDir = new File(datadir.resolve("interpro").toString())
        if (!_interproDir.exists() || !_interproDir.isDirectory()) {
            error = "No 'interpro' directory was found in the data directory ${datadir}\n" +
                    "Please ensure that the data dir is correctly populated or use --download"
            return [_interproRelease, error]
        }

        if (interpro == "latest") {
            // Use the latest interpro release in datadir/interpro
            def _dirs = _interproDir.listFiles()
                    .findAll { it.isDirectory() }
                    .collect { it.name.toFloat() } // what if the user put a dir in that can't be converted to a float
            if (_dirs.isEmpty()) {
                error = "No InterPro release directories were found in ${datadir}/interpro" +
                        "Please ensure that the data dir is correctly populated or use --download"
                return [_interproRelease, error]
            }
            _interproRelease = _dirs.max()
        } else {
            // Use the user specified InterPro release
            _interproRelease = interpro
            def _interproReleaseDir = new File(datadir.resolve("interpro/${_interproRelease}").toString())
            if (!_interproReleaseDir.exists() || !_interproReleaseDir.isDirectory()) {
                error = "Could not find the InterPro release dir at ${_interproReleaseDir}\n" +
                        "Please check the correct InterPro release was provided and that the data dir is correctly populated or use --download"
                return [_interproRelease, error]
            }
        }
        return [_interproRelease.toString(), error]
    }

    static String getDatabaseVersion(String database, String databaseJson) {
        // Retrieve the database version from xrefs/databases.json
        // var database may be a Path or string
        JsonSlurper jsonSlurper = new JsonSlurper()
        File file = new File(databaseJson)
        def memberDbReleases = jsonSlurper.parse(file)
        return memberDbReleases[database]
    }

    static httpRequest(String urlString, String data, int maxRetries, boolean verbose, log) {
        int attempts = 0
        boolean isPost = data != null && data.length() > 0
        HttpURLConnection connection = null
        while (attempts < maxRetries + 1) {
            attempts++

            try {
                URL url = new URL(urlString)
                connection = (HttpURLConnection) url.openConnection()
                connection.connectTimeout = 5000
                connection.readTimeout = 5000
                if (isPost) {
                    connection.doOutput = true
                    connection.requestMethod = "POST"
                    connection.setRequestProperty("Content-Type", "application/json")
                    connection.setRequestProperty("Accept", "application/json")
                    connection.outputStream.withWriter("UTF-8") { writer ->
                        writer << data
                    }
                    connection.outputStream.flush()
                } else {
                    connection.requestMethod = "GET"
                }

                int responseCode = connection.responseCode
                if (responseCode >= 200 && responseCode < 300) {
                    String responseText = connection.inputStream.getText("UTF-8")
                    return new JsonSlurper().parseText(responseText)
                } else {
                    if (verbose) {
                        def errorMsg = connection.errorStream ? connection.errorStream.getText("UTF-8") : ""
                        log.warn("Received HTTP ${responseCode} for ${urlString}")
                    }
                    return null
                }
            } catch (java.net.ConnectException | java.net.UnknownHostException e) {
                if (verbose) {
                    log.warn("Connection error for ${urlString}: ${e.message}")
                }
                return null
            } catch (java.net.SocketTimeoutException e) {
                if (verbose) {
                    log.warn("Timeout error for ${urlString}; attempt ${attempts} / ${maxRetries + 1}")
                }
            } catch (Exception e) {
                if (verbose) {
                    log.warn("Unexpected error for ${urlString}: ${e.message}")
                }
                return null
            } finally {
                connection?.disconnect()
            }
        }
    }

    static String sanitizeURL(String url) {
        return url.replaceAll(/\/+$/, '')
    }

    static List<String> getFtpMd5(String database, String dbRelease) {
        // Retrieve the precompute md5 hash from the EBI ftp
        String error = ""
        String md5Hash = ""
        String url = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/${database}/${database}-${dbRelease}.tar.gz.md5"
        try {
            md5Hash = new URL(url).text.replace("${database}-{$dbRelease}.tar.gz", "").trim()
        } catch (IOException e) {
            error = "Error: Unable to retrieved MD5 hash from the EBI FTP at ${url}"
        }
        return [md5Hash, error]
    }
}