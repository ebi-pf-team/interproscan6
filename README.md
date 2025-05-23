# InterProScan 6

[![Nextflow](https://img.shields.io/badge/%E2%89%A524.10.4-0dc09d?style=flat&logo=nextflow&label=nextflow&labelColor=f5fafe)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/docker-1D63ED?logo=docker&logoColor=1D63ED&label=run%20with&labelColor=f5fafe)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/singularity-1E4383?label=run%20with&labelColor=f5fafe)](https://sylabs.io/docs/)

[InterPro](http://www.ebi.ac.uk/interpro/) is a database that brings together predictive information on protein function from multiple partner resources. It provides an integrated view of the families, domains and functional sites to which a given protein belongs.

**InterProScan** is the command‑line tool that allows you to scan protein or nucleotide sequences against the InterPro member‑database signatures in a single workflow. Researchers with novel sequences can use InterProScan to annotate their data with family classifications, domain architectures and site predictions.

> [!CAUTION]
> InterProScan 6 is under active development and may be unstable.

## Installation

Before you begin, you will need:

* [Nextflow](https://www.nextflow.io/) 24.10.04 or later
* A container runtime. The currently supported are:
    * [Docker](https://www.docker.com/)
    * [SingularityCE](https://sylabs.io/singularity/)
    * [Apptainer](https://apptainer.org/)

> [!IMPORTANT]  
> Phobius, SignalP and DeepTMHMM also require separate licenses and data downloads. See [Licensed analyses](#licensed-analyses) below.

## Usage

### Running the built‑in test

To test InterProScan, run the following command:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile test,docker \
  --datadir data \
  --interpro latest \
  --download
```

Explanation of parameters:

* `-r 6.0.0-alpha`: Specify the version of InterProScan to use, in this case version `6.0.0-alpha`. For consistent and reproducible results, we strongly recommend always specifying a version.
* `-profile test,docker`:
  * `test`: Run InterProScan with a small example FASTA file included in the workflow.
  * `docker`: Execute tasks in Docker containers.
* `--datadir data`: Set `data` as the directory for storing all required InterPro and member database files. The directory is created automatically if it doesn't exist.
* `--interpro latest`: Use the most recent InterPro release.
* `--download`: Download any missing metadata and database files into the specified `--datadir`

> [!TIP]
> To use a specific version of InterPro, use the release number, e.g. `--interpro 105.0`.

After completion, you'll find three output files in your working directory:

* `test.faa.json`: full annotations (JSON)
* `test.faa.tsv`: tabular summary (TSV)
* `test.faa.xml`: full annotations (XML)

The JSON and XML outputs are more comprehensive; the TSV is a concise summary.

### Scanning your own sequences

To annotate your own sequences FASTA file, omit the `test` profile and specify `--input`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.fasta
```

For nucleotide sequences, add `--nucleic`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.fna \
  --nucleic
```

### Selecting specific analyses

To run only certain analyses (for example, Pfam and MobiDB‑lite):

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile test,docker \
  --datadir data \
  --applications Pfam,MobiDB-lite
```

Analyses are case‑insensitive and ignore hyphens, so both `MobiDB-lite` and `mobidblite` work.

The available analyses are:

* `AntiFam`
* `CATH-Gene3D`
* `CATH-FunFam`
* `CDD`
* `COILS`
* `DeepTMHMM`
* `HAMAP`
* `MobiDB-lite`
* `NCBIFAM`
* `PANTHER`
* `Phobius`
* `Pfam`
* `PIRSF`
* `PIRSR`
* `PRINTS`
* `PROSITE-patterns`
* `PROSITE-profiles`
* `SFLD`
* `SMART`
* `SUPERFAMILY`
* `SignalP-Euk`
* `SignalP-Prok`

> [!TIP]
> If you only run COILS, DeepTMHMM, MobiDB‑lite, Phobius, SignalP‑Euk or SignalP‑Prok, you don't need `--datadir`, `--interpro`, or `--download`.

### Running InterProScan on an HPC cluster with Slurm

To run InterProScan on your institute's Slurm cluster, use the `slurm` profile. This ensures that each task in the pipeline is submitted as a job to the Slurm scheduler.

Most HPC systems do not support Docker, but they often support Singularity or Apptainer for containerized execution. Include the appropriate profile (`singularity` or `apptainer`).

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile singularity,slurm,test \
  --datadir /path/to/data
```

> [!IMPORTANT]
> The directory specified by `--datadir` must be accessible from all cluster nodes. This usually means it should be located on a shared network file system (e.g. NFS or Lustre).

## Licensed analyses

DeepTMHMM, Phobius and SignalP contain licensed components and are disabled by default. To enable and execute these analyses, you first need to obtain their license, and data.

### Obtaining licensed components

For each analyses, you need to request a license, then download and extract an archive that contains the data and machine learning model(s) required to run the analysis.

#### DeepTMHMM 1.0

Request a standalone copy of DeepTMHMM 1.0 by sending an email to <licensing@biolib.com>, then extract it when you receive it:

```sh
unzip -q DeepTMHMM-v1.0.zip
```

and make note of the full path to the package:

```sh
echo "${PWD}/DeepTMHMM
```

#### Phobius 1.01

Download a copy of Phobius 1.01 [from Erik Sonnhammer's website](https://software.sbc.su.se/phobius.html), then extract it:

```sh
tar -zxf phobius101_linux.tgz
```

and make note of the full path to the package:

```sh
echo "${PWD}/phobius"
```

#### SignalP 6.0

SignalP 6.0 supports two modes: a slow one that uses the full model, and a fast one that uses a distilled version of the full mode. InterProScan suppports both, but only one at a time. The fast mode is recommended for most users.

You need a license for each of these models:

* download the distilled model: [fast](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast)
* download the full model: [slow_sequential](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=slow_sequential)

Extract the archive:

```sh
tar -zxf signalp-6.0h.fast.tar.gz
```

and make note of the full path to the package:

```sh
echo "${PWD}/signal6p_fast"
```

### Executing licensed analyses

Create a Nextflow config (e.g. `licensed.conf`) that points to each tool directory:

```groovy
params {
    appsConfig {
        deeptmhmm {
            dir = "/full/path/to/DeepTMHMM"
        }
        phobius {
            dir = "/full/path/to/phobius"
        }
        signalp_euk {
            dir = "/full/path/to/signal6p_fast"
        }
        signalp_prok {
            dir = "/full/path/to/signal6p_fast"
        }
    }
}
```

And run with your config:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-alpha \
  -profile docker \
  -c licensed.conf \
  --input /path/to/sequences.fasta \
  --applications deeptmhmm,phobius,signalp_euk,signalp_prok \
  --offline
```

> [!WARNING]  
> As DeepTMHMM 1.0 and SignalP 6.0 predictions are not yet available in the [Matches API](https://www.ebi.ac.uk/interpro/matches/api/), the pre-calculated matches lookup must be disabled with `--offline`.

> [!WARNING]  
> Phobius does not support certain non-standard or ambiguous residues. Any sequence containing pyrrolysine (one-letter code `O`), Asx (Asp/Asn ambiguity, `B`), Glx (Glu/Gln ambiguity, `Z`) or Xle (Leu/Ile ambiguity, `J`) will be skipped by Phobius but will continue to be processed normally by all other applications.

> [!NOTE]  
> Running both `signalp_euk` and `signalp_prok` will execute SignalP twice, once with eukaryotic post-processing and once without. Choose the mode best suited to your dataset.

## Documentation

Our full documentation is available on [ReadTheDocs](https://interproscan-docs.readthedocs.io/en/v6/).

## Citation

If you use InterPro in your work, please cite the following publication:

> Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A. **InterPro: the protein sequence classification resource in 2025**. *Nucleic Acids Res*. 2025 Jan;53(D1):D444-D456. [doi: 10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).

## Support

For further assistance, please [create an issue](https://github.com/ebi-pf-team/interproscan6/issues) or [contact us](http://www.ebi.ac.uk/support/interproscan).
