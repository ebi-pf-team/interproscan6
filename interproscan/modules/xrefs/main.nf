process GOTERMS {
    label 'xref'

    input:
    val ipr2go
    val go_info
    val matchesinterpro

    output:
    val matches2go

    script:
    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    def matches_info = matchesinterpro[1]
    def match_key2go = [:]
    matches_info.each { match_key, interpro_key ->
        try {
            def go_ids = ipr2go[interpro_key]
            for (String go_id : go_ids) {
                def go_dict = [
                    "name": go_info[go_id][1],
                    "databaseName": GO_PATTERN[go_info[go_id][0]],
                    "id": go_id
                ]
                if (!match_key2go.containsKey(match_key)) {
                    match_key2go[match_key] = []
                }
                match_key2go[match_key] << go_dict
            }
        } catch (Exception e) {
            match_key2go[match_key] = []
        }
    }
    return match_key2go
}

process PATHWAYS {
    label 'xref'

    input:
    val ipr2pa
    val pa_info
    val matchesinterpro

    output:
    val matches2pa

    script:
    def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
    ]

    def matches_info = matchesinterpro[1]
    def match_key2pa = [:]
    matches_info.each { match_key, interpro_key ->
        try {
            def go_ids = ipr2pa[interpro_key]
            for (String go_id : pa_ids) {
                def pa_dict = [
                    "name": pa_info[pa_id][1],
                    "databaseName": GO_PATTERN[pa_info[pa_id][0]],
                    "id": pa_id
                ]
                if (!match_key2pa.containsKey(match_key)) {
                    match_key2pa[match_key] = []
                }
                match_key2pa[match_key] << pa_dict
            }
        } catch (Exception e) {
            match_key2pa[match_key] = []
        }
    }
    return match_key2pa
}
