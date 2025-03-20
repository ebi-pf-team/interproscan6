// Class and methods for interacting with InterPro servers and data

import groovy.json.JsonSlurper
import java.net.URL
import com.fasterxml.jackson.databind.ObjectMapper

class InterPro {
    static String getDatabaseVersion(database, directory) {
        // Retrieve the database version from xrefs/databases.json
        // database may be a Path or string
        // TODO: change to JsonSlurper when we migrate to databases.json
        ObjectMapper objectMapper = new ObjectMapper();
        File file = new File(new File(directory.toString(), "xrefs"), "entries.json")
        Map<String, Object> metadata = objectMapper.readValue(file, Map.class);
        return metadata.databases[database]
    }

    static Map getInfo(String baseUrl, log) {
        return httpRequest("${sanitizeURL(baseUrl)}/info", null, 0, true, log)
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

    static String getMatchesApiUrl(matchesApiUrl, lookupServiceUrl, datadir, workflowManifest, log) {
        String _matchesApiUrl = matchesApiUrl ?: lookupServiceUrl
        Map info = getInfo(_matchesApiUrl, log)

        if (info == null) {
            log.warn("An error occurred while querying the Matches API; analyses will be run locally")
            _matchesApiUrl = null
        } else {
            def apiVersion = info.api ?: "X.Y.Z"
            def majorVersion = apiVersion.split("\\.")[0]

            if (majorVersion != "0") {
                log.warn("${workflowManifest.name} ${workflowManifest.version}" +
                        " is not compatible with the Matches API at ${_matchesApiUrl};" +
                        " analyses will be run locally")
                _matchesApiUrl = null
            } else if (datadir) {
                def iprVersion = getDatabaseVersion("InterPro", datadir)
                if (iprVersion != info.release) {
                    log.warn("InterPro version mismatch (local: ${iprVersion}, Matches API: ${info.release});" +
                            " pre-calulated matches will not be retrieved, and analyses will be run locally")
                    _matchesApiUrl = null
                }
            }
        }

        return _matchesApiUrl
    }

    static String sanitizeURL(String url) {
        return url.replaceAll(/\/+$/, '')
    }
}