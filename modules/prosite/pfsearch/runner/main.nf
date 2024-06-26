process PFSEARCH_RUNNER {
    /*
    We us an inhouse python script to coordinate running pfsearch for 
    all provided PROSTITE Profiles
    */
    container 'docker.io/sibswiss/pftools'
    label 'prosite_pfsearch_runner'

    input:
        tuple path(fasta), path(models_dir), val(release), val(switches)

    output:
        path "${release}_prosite_profiles.out"

    script:
    """
    $projectDir/scripts/prosite/run_prosite.py \
        ${models_dir} \
        ${fasta} \
        ${release}_prosite_profiles.out \
        pfsearchV3 \
        ${swithces} 
    """
}