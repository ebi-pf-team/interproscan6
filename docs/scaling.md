# Scaling InterProScan

This page covers the built-in profiles for scaling InterProScan 6 on HPC schedulers and cloud batch systems:

- HPC schedulers: `slurm`, `lsf`
- cloud batch systems: `googlebatch`, `awsbatch`

For scalable runs, InterProScan expects a container runtime that matches the platform:

- on HPC systems, this will usually be `singularity` or `apptainer`
- on Google Cloud Batch and AWS Batch, use Docker-based execution

## HPC clusters

InterProScan provides built-in scheduler profiles for:

- `slurm`
- `lsf`

These profiles submit tasks to the scheduler instead of running everything on the local machine.

Recommended profile combinations on clusters include:

!!! warning

    The directory specified by `--datadir` and the Nextflow working directory must be accessible from compute nodes. In practice, this usually means using a shared filesystem such as NFS, Lustre, or GPFS.

Example commands:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile singularity,slurm \
  --input /path/to/sequences.faa \
  --datadir /shared/interproscan/data
```

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile apptainer,lsf \
  --input /path/to/sequences.faa \
  --datadir /shared/interproscan/data
```

## Cloud platforms

InterProScan provides built-in cloud batch profiles for:

- `awsbatch`: AWS Batch
- `googlebatch`: Google Batch

These profiles define the executor behavior inside the workflow, but you still need to provide site-specific cloud configuration at runtime with an additional Nextflow config file passed by `-c`.

### AWS Batch

Use the `awsbatch` profile and pass your AWS-specific config on the command line. Cloud batch runs should use Docker-based execution.

Example `aws.config`:

```groovy
process {
    queue       = '<aws-batch-queue>'
}

aws {
    accessKey   = '<access-key>'
    secretKey   = '<secret-key>'
    region      = '<aws-region>'

    batch {
    	maxSpotAttempts = 10 // (1)!
    }
}
```

1. Set the maximum number of retries to 10, when using Spot instances

!!! tip

    There are multiple ways to provide security credentials and the region at runtime. See [AWS security credentials](https://docs.seqera.io/nextflow/aws#aws-security-credentials) on the Nextflow documentation.

Example command:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker,awsbatch \
  -c aws.config \
  -bucket-dir s3://<bucket>/interproscan/work \
  --input /path/to/sequences.faa \
  --datadir s3://<bucket>/interproscan/data
```

### Google Cloud Batch

Use the `googlebatch` profile and pass your GCP-specific config on the command line. Cloud batch runs should use Docker-based execution.

Example `gcp.config`:

```groovy
google {
  project = '<your-gcp-project>'
  location = '<gcp-region>'

  batch {
  	spot = true// (1)!
    maxSpotAttempts = 10 // (2)!
  }
}
```

1. Enable the use of virtual machines
2. Set the maximum number of retries to 10, when using `spot` is set to `true`

Example command:

```bash
export GOOGLE_APPLICATION_CREDENTIALS="/path/to/creds.json"
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -profile docker,googlebatch \
  -c gcp.config \
  -bucket-dir gs://<bucket>/interproscan/work \
  --input /path/to/sequences.faa \
  --datadir gs://<bucket>/interproscan/data
```

## GPU-enabled runs

Use `--use-gpu` to enable GPU acceleration for compatible deep-learning-based analyses that have been selected and configured.

The workflow automatically handles the basic GPU request for supported execution environments:

- Docker: adds `--gpus all`
- Singularity and Apptainer: add `--nv`
- Slurm: requests `--gres=gpu:1`
- LSF: requests `-gpu "num=1"`
- AWS Batch: requests one accelerator
- Google Batch: requests one accelerator

You are still responsible for choosing the GPU type where the platform requires it.

Typical scheduler overrides look like this for Slurm:

```groovy
process {
  withLabel: use_gpu {
    clusterOptions = '--gres=gpu:a100:1'
  }
}
```

And for LSF:

```groovy
process {
  withLabel: use_gpu {
    clusterOptions = '-gpu "num=1:gmodel=A100"'
  }
}
```

Add such overrides in an extra Nextflow config file and pass it with `-c your_cluster.config`.

Example command:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.1 \
  -c your_cluster.config \
  -profile singularity,slurm \
  --input /path/to/sequences.faa \
  --datadir /shared/interproscan/data \
  --run-ml \
  --use-gpu
```
