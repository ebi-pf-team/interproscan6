process PFSEARCH_RUNNER {
    /*
    We use an inhouse python script to coordinate running pfsearch for 
    all provided PROSTITE Profiles
    */
    // container 'docker.io/sibswiss/pftools'
    label 'prosite_pfsearch_runner'

    input:
        tuple path(fasta), path(models_dir), val(release), val(switches), path(blacklist_file)

    output:
        path "${release}._.prosite_profiles.out"
        path blacklist_file

    // change to use just pfsearchV3 not the full project dir path
    // when we move over to the single docker file
    script:
    """
    python3 $projectDir/scripts/members/prosite/run_pfsearchv3.py \
        ${models_dir} \
        ${fasta} \
        ${release}._.prosite_profiles.out \
        "/opt/pftools/var/lib/pftools/bin/pfsearchV3" \
        ${switches}
    """
}