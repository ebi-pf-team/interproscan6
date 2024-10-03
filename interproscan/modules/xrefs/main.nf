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
    groovy -cp . XRefs.groovy XRefs match_key2go '${ipr2go}' '${go_info}' '${matchesinterpro}' > matches2go
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
    groovy -cp . XRefs.groovy XRefs match_key2pa '${ipr2pa}' '${pa_info}' '${matchesinterpro}' > matches2pa
    """
}
