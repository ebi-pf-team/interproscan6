process PFSEARCH_RUNNER {
    /*
    We use an inhouse python script to coordinate running pfsearch for
    all provided PROSTITE Profiles
    */
    label 'prosite_pfsearch_runner'

    input:
        tuple val(meta), path(fasta)
        path models_dir
        path blacklist_file

    output:
        path "prosite_profiles.out"
        path blacklist_file

    script:
    """
    touch prosite_profiles.out
    python3 $projectDir/interproscan/scripts/members/prosite/run_pfsearchv3.py \
        ${models_dir} \
        ${fasta} \
        prosite_profiles.out \
        "/opt/pftools/pfsearchV3" \
        -f -o 7 -t 4
    """
}


process PFSEARCH_PARSER {
    label 'analysis_parser'

    input:
        tuple val(meta), path(pfsearch_out)
        path blacklist_file

    output:
        tuple val(meta), path("${pfsearch_out}-filtered.json")

    exec:
//     pfsearch_parser.py \
//         ${pfsearch_out} \
//         ${pfsearch_out}-filtered.json \
//         ${blacklist_file}
}
