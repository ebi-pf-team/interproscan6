import java.net.URI
import groovy.json.JsonSlurper

class HTTPRequest {
    static fetch(String urlString, String data, int maxRetries, log) {
        int attempts = 0
        boolean isPost = data != null && data.length() > 0
        HttpURLConnection connection = null
        while (attempts < maxRetries + 1) {
            attempts++

            try {
                URL url = new URI(urlString).toURL()
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
                    if (log != null) {
                        def errorMsg = connection.errorStream ? connection.errorStream.getText("UTF-8") : ""
                        log.warn("Received HTTP ${responseCode} for ${urlString}")
                    }
                    return null
                }
            } catch (java.net.ConnectException | java.net.UnknownHostException e) {
                if (log != null) {
                    log.warn("Connection error for ${urlString}: ${e.message}")
                }
                return null
            } catch (java.net.SocketTimeoutException e) {
                if (log != null) {
                    log.warn("Timeout error for ${urlString}; attempt ${attempts} / ${maxRetries + 1}")
                }
            } catch (Exception e) {
                if (log != null) {
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