process GOTERMS {
    label 'xref'

    input:
    val ipr2go
    val go_info
    val matches2entries

    output:
    val matches2go

    script:
    """
    def GO_PATTERN = [
        "P": "BIOLOGICAL_PROCESS",
        "C": "CELLULAR_COMPONENT",
        "F": "MOLECULAR_FUNCTION"
    ]

    matches2go = matches_info.collectEntries { seq_id, match_info ->
        match_info.each { match_key, data ->
            if (data["entry"]) {
                interpro_key = data["entry"]["accession"]
                try {
                    def go_ids = ipr2go[interpro_key]
                    for (String go_id : go_ids) {
                        def go_dict = [
                            "name": go_info[go_id][0],
                            "databaseName": "GO",
                            "category": GO_PATTERN[go_info[go_id][1]],
                            "id": go_id
                        ]
                        match_info[match_key]["entry"]["goXRefs"] << go_dict
                    }
                } catch (Exception e) {
                    // pass (no GO terms for this interpro_key)
                }
            }
        }
        return match_info
    }
    """
}

process PATHWAYS {
    label 'xref'

    input:
    val ipr2pa
    val pa_info
    val matches2entries

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

    matches2pa = matches_info.collectEntries { seq_id, match_info ->
        match_info.each { match_key, data ->
            if (data["entry"]) {
                interpro_key = data["entry"]["accession"]
                try {
                    pa_ids = ipr2pa[interpro_key]
                    for (String pa_id : pa_ids) {
                        pa_dict = [
                            "name": pa_info[pa_id][1],
                            "databaseName": PA_PATTERN[pa_info[pa_id][0]],
                            "id": pa_id
                        ]
                        match_info[match_key]["entry"]["pathwayXRefs"] << pa_dict
                    }
                } catch (Exception e) {
                    // pass (no Pathways for this interpro_key)
                }
            }
        }
        return matches_info
    }
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
