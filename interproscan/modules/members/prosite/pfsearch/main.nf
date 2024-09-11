// Scripts for Prosite Profiles

process PFSEARCH_RUNNER {
    /*
    We use an inhouse python script to coordinate running pfsearch for
    all provided PROSTITE Profiles
    */
    label 'prosite_pfsearch_runner'

    input:
        tuple path(fasta), path(models_dir), val(switches), path(blacklist_file)

    output:
        path "prosite_profiles.out"
        path blacklist_file

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prosite/run_pfsearchv3.py \
        ${models_dir} \
        ${fasta} \
        prosite_profiles.out \
        "/opt/pftools/pfsearchV3" \
        ${switches}
    """
}


process PFSEARCH_PARSER {
    label 'analysis_parser'

    input:
        path pfsearch_out
        path blacklist_file

    output:
        path "${pfsearch_out}-filtered.json"

    script:
    """
    python3 $projectDir/interproscan/scripts/members/prosite/pfsearch_parser.py \
        ${pfsearch_out} \
        ${pfsearch_out}-filtered.json \
        ${blacklist_file}
    """
}
