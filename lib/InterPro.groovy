// Class and methods for interacting with InterPro servers and data

import groovy.json.JsonSlurper
import java.net.URL
import com.fasterxml.jackson.databind.ObjectMapper
import InterProScan

class InterPro {
    static final def FTP_URL = "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6"

    static boolean checkCompatibility(Float iprscan, String interpro) {
        // Check InterPro and InterProScan versions are compatible
        def noConnWarning = ""
        def latestReleaseWarning = ""
        def error = ""
        def releasesURL = InterPro.sanitizeURL("${FTP_URL}/${iprscanURL}/releases.json")
        Map compatibleReleases = InterPro.httpRequest(releasesURL, null, 0, true, log)
        if (compatibleReleases == null) {
            noConnWarning = "An error occurred while querying the EBI ftp to assess the compatibility of\n"+
                     "InterPro release '${interpro}' with InterProScan version '${iprscan}'\n"+
                     "Proceeding with using InterPro release '${interpro}'. Warning InterProScan may not behave as expected"
        } else if (!compatibleReleases[iprscan].contains(interpro)) {
            error = "InterPro release '$interpro' is not compatiable with "
        } else if (interpro < compatibleReleases[iprscan].max()) {
            latestReleaseWarning = "A later compatible InterPro release '${compatibleReleases[iprscan].max()}' is availble.\n"+
                                   "Tip: Use the '--download' option to automatically fetch the latest compatible InterPro release"
        }
        return [noConnWarning, latestReleaseWarning, error]
    }

    static getInterproRelease(Path datadir, String interpro) {
        // Checks the InterProDir exists and then gets the specified interpro release number
        def _interproRelease = null
        def _interproDir = ""
        def error = ""

        (_interproDir, error) = InterProScan.validateInterproDir(datadir)
        if (error) {  // no datadir/interpro dir found
            return [_interproRelease, error]
        }

        if (interpro == "latest") {
            // Use the latest interpro release in datadir/interpro
            def _dirs = _interproDir.listFiles()
                    .findAll { it.isDirectory() }
                    .collect { it.name.toFloat } // what if the user put a dir in that can't be converted to a float
            if (_dirs.isEmpty()) {
                error = "No InterPro release directories were found in ${datadir}/interpro" +
                        "Please ensure that the data dir is correctly populated or use --download"
                return [_interproRelease, error]
            }
            _interproRelease = _dirs.max()
        } else {
            // Use the user specified InterPro release
            _interproRelease = interpro
        }
        return [_interproRelease, error]
    }

    static String getDatabaseVersion(database, directory) {
        // Retrieve the database version from xrefs/databases.json
        // database may be a Path or string
        // TODO: change to JsonSlurper when we migrate to databases.json
        ObjectMapper objectMapper = new ObjectMapper();
        File file = new File(new File(directory.toString(), "xrefs"), "entries.json")
        Map<String, Object> metadata = objectMapper.readValue(file, Map.class);
        return metadata.databases[database]
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
}