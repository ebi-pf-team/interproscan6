process GOTERMS {
    label 'xref'

    input:
    val ipr2go
    val go_info
    val matchesinterpro

    output:
    val matches2go

    script:
    """
    import groovy.json.JsonOutput

    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    def matches2go = [:]
    matchesinterpro[1].each { match_key, interpro_key ->
        try {
            go_ids = ipr2go[interpro_key]
            for (String go_id : go_ids) {
                go_dict = [
                    "name": go_info[go_id][1],
                    "databaseName": GO_PATTERN[go_info[go_id][0]],
                    "id": go_id
                ]
                matches2go[match_key] << go_dict
            }
        } catch (Exception e) {
            // pass (no GoTerms for this interpro_key)
        }
    }
    return matches2go
    """
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
    """
    def PA_PATTERN = [
    "t": "MetaCyc",
    "w": "UniPathway",
    "k": "KEGG",
    "r": "Reactome"
    ]

    def matches2pa = [:]
    matchesinterpro[1].each { match_key, interpro_key ->
        try {
            pa_ids = ipr2pa[interpro_key]
            for (String pa_id : pa_ids) {
                pa_dict = [
                    "name": pa_info[pa_id][1],
                    "databaseName": PA_PATTERN[pa_info[pa_id][0]],
                    "id": pa_id
                ]
                matches2pa[match_key] << pa_dict
            }
        } catch (Exception e) {
            // pass (no Pathways for this interpro_key)
        }
    }
    return matches2pa
    """
}

process PANTHER_XREFS {
    label 'xref'

    input:
    val paint_annotations
    val matches2entries

    output:
    val matches2pthrxrefs

    script:
    """
    matches2pthrxrefs = matches_info.collectEntries { seq_id, match_info ->
        match_info.each { match_key, data ->
            def acc_subfamily = data["accession"]
            try {
                match_info[match_key]["entry"]["subfamily_name"] = entries_info[acc_subfamily]["name"]
                match_info[match_key]["entry"]["subfamily_description"] = entries_info[acc_subfamily]["description"]
                match_info[match_key]["entry"]["subfamily_type"] = entries_info[acc_subfamily]["type"]
            } catch (Exception e) {
                log.info "No subfamily for ${acc_subfamily}"
            }
            node_data = paint_annotations[data["node_id"]]
            match_info[seq_id][sig_acc]["proteinClass"] = node_data[2]
            match_info[seq_id][sig_acc]["graftPoint"] = node_data[3]
        }
        return match_info
    }

    """
}
