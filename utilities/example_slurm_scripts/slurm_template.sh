#!/bin/bash
#SBATCH --job-name=ips6-test
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=30:00
#SBATCH --mem=20G

cd $SLURM_SUBMIT_DIR

# add code to load container runtime and nextflow software

nextflow run interproscan.nf \
    -profile slurm,singularity \
    -work-dir /lscratch \
    --input utilities/test_files/best_to_test.fasta \
    --disable_precalc \
    --applications antifam,ncbifam,pfam
