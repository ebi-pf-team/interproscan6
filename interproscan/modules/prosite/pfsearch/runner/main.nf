process PFSEARCH_RUNNER {
    /*
    We use an inhouse python script to coordinate running pfsearch for
    all provided PROSTITE Profiles
    */
    // container 'docker.io/sibswiss/pftools'
    label 'prosite_pfsearch_runner'

    input:
        tuple path(fasta), path(models_dir), val(switches), path(blacklist_file)

    output:
        path "prosite_profiles.out"
        path blacklist_file

    // change to use just pfsearchV3 not the full project dir path
    // when we move over to the single docker file
    script:
    """
    touch prosite_profiles.out
    python3 $projectDir/interproscan/scripts/members/prosite/run_pfsearchv3.py \
        ${models_dir} \
        ${fasta} \
        prosite_profiles.out \
        "/opt/pftools/pfsearchV3" \
        ${switches}
    """
}
