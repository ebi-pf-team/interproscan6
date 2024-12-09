process RUN_TMHMM {
    label 'medium', 'tmhmm_runner'
    stageInMode 'copy'

    input:
    tuple val(meta), path(fasta)
    path tmhmm_dir

    output:
    tuple val(meta), path("outdir")

    script:
    """
    # deeptmhmm has a hard coded assumption it is being run within its dir
    cd ${tmhmm_dir}
    python3 predict.py \
        --fasta ../${fasta} \
        --output-dir ../outdir
    """
}
