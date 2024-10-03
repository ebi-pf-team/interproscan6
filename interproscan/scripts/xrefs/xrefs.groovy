class XRefs {
    static final def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]
    static final def PA_PATTERN = [
        "t": "MetaCyc",
        "w": "UniPathway",
        "k": "KEGG",
        "r": "Reactome"
    ]

    static Map match_key2go(Map ipr2go, Map go_info, Map matchesinterpro) {
        def match_key2go = [:]
        matchesinterpro[1].each { match_key, interpro_key ->
            try {
                go_ids = ipr2go[interpro_key]
                for (String go_id : go_ids) {
                    def go_dict = [
                        "name": go_info[go_id][1],
                        "databaseName": GO_PATTERN[go_info[go_id][0]],
                        "id": go_id
                    ]
                    match_key2go[match_key] = match_key2go.get(match_key, [])
                    match_key2go[match_key] << go_dict
                }
            } catch (Exception e) {
                match_key2go[match_key] = []
            }
        }
        return match_key2go
    }

    static Map match_key2pa(Map ipr2pa, Map pa_info, Map matchesinterpro) {
        def match_key2pa = [:]
        matchesinterpro[1].each { match_key, interpro_key ->
            try {
                def pa_ids = ipr2pa[interpro_key]
                for (String pa_id : pa_ids) {
                    pa_dict = [
                        "name": pa_info[pa_id][1],
                        "databaseName": PA_PATTERN[pa_info[pa_id][0]],
                        "id": pa_id
                    ]
                    match_key2pa[match_key] = match_key2pa.get(match_key, [])
                    match_key2pa[match_key] << pa_dict
                }
            } catch (Exception e) {
                match_key2pa[match_key] = []
            }
        }
        return match_key2pa
    }
}
