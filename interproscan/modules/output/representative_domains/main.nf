process REPRESENTATIVE_DOMAINS {
    label 'io'

    input:
    path matches

    output:
    path "repr_domains_${matches}"

    script:
    """
    python3 $projectDir/interproscan/scripts/output/representative_domains.py ${matches} > repr_domains_${matches}
    """
}
