import uk.ac.ebi.interpro.MatchesApiClient

workflow PREPARE_APPLICATIONS {
    take:
    apps             // list[str], names of applications to prepare
    no_matches_api   // boolean, whether to retrieve matches from the Matches API
    matches_api_url  // str, url to the Matches API
    name             // str, name of the workflow
    version          // str, version of the workflow

    main:
    (matches_api_apps, local_only_apps, api_version, error) = MatchesApiClient.prepareLookup(
        apps,
        no_matches_api,
        matches_api_url,
        name,
        version
    )
    if (error) {
        log.warn error
    } else if (!no_matches_api && !local_only_apps.isEmpty()) {
        log.warn "The following analyses are not available in the Matches API: " +
                local_only_apps.join(", ") +
                ". They will be executed locally."
    }

    emit:
    local_only_apps  // list: list of application to that are not in the matches API
    matches_api_apps // list: list of applications that are in the matches API
    api_version      // str: version of the matches API
}