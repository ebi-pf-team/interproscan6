class Lookup {
    static prepareLookup(List<String> apps, Boolean noMatchesApi, String matchesApiUrl, workflowManifest) {
        List<String> allMatchesApiApps = []
        List<String> matchesApiApps = []
        List<String> localOnlyApps  = []
        String apiInterproVersion = null
        String apiVersion
        String error

        if (!noMatchesApi && matchesApiUrl) { // url check for unit tests
            Map info = HTTPRequest.fetch("${HTTPRequest.sanitizeURL(matchesApiUrl)}/info".toString(), null, 0, true)
            if (info == null) {
                error = "An error occurred while querying the Matches API;" +
                        " analyses will be run locally"
            } else {
                apiVersion = info.api ?: "X.Y.Z"
                def majorVersion = apiVersion.split("\\.")[0]
                if (majorVersion != "0") {
                    error = "${workflowManifest.name} ${workflowManifest.version}" +
                            " is not compatible with the Matches API at ${matchesApiUrl};" +
                            " analyses will be run locally"
                } else {
                    apiInterproVersion = info.release
                    if (info.analyses) {
                        allMatchesApiApps.addAll(info.analyses*.name.collect { it.toLowerCase().replaceAll("[-\\s]", "") })
                        apps.each { app ->
                            (allMatchesApiApps.contains(app.toLowerCase().replaceAll("[-\\s]", "")) ? matchesApiApps : localOnlyApps) << app
                        }
                    } else {
                        error = "Could not retrieve the list of available analyses from the Matches API; " +
                                " analyses will be run locally"
                    }
                }
            }   
        }

        if (error || matchesApiApps.isEmpty() || noMatchesApi) {
            localOnlyApps  = apps.clone() as List<String>
        }

        return [matchesApiApps, localOnlyApps, apiInterproVersion, error]
    }
}
