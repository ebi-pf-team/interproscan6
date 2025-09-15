# InterProScan 6

[![Nextflow](https://img.shields.io/badge/%E2%89%A524.10.4-0dc09d?style=flat&logo=nextflow&label=nextflow&labelColor=f5fafe)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/docker-1D63ED?logo=docker&logoColor=1D63ED&label=run%20with&labelColor=f5fafe)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/singularity-1E4383?label=run%20with&labelColor=f5fafe)](https://sylabs.io/docs/)

[InterPro](http://www.ebi.ac.uk/interpro/) is a database that brings together predictive information on protein function from multiple partner resources. It provides an integrated view of the families, domains and functional sites to which a given protein belongs.

**InterProScan** is the command‑line tool that allows you to scan protein or nucleotide sequences against the InterPro member‑database signatures in a single workflow. Researchers with novel sequences can use InterProScan to annotate their data with family classifications, domain architectures and site predictions.

> [!WARNING]
> **InterProScan 6 is now in beta**.
> The core functionality is complete, and results have been benchmarked against InterProScan 5 with high similarity across all analyses.
>
> We encourage users to test this version and report any issues or suggestions. Minor changes and refinements may still occur before the final release.

## Installation

Before you begin, install:

* [Nextflow](https://www.nextflow.io/) 25.04 or later
* A container runtime. The currently supported are:
    * [Docker](https://www.docker.com/)
    * [SingularityCE](https://sylabs.io/singularity/)
    * [Apptainer](https://apptainer.org/)

You don't need anything else, Nextflow will download the workflow from GitHub, 
and required data are automatically downloaded when running InterProScan.

> [!IMPORTANT]  
> Phobius, SignalP and DeepTMHMM require separate licenses and downloads. See [Licensed analyses](#licensed-analyses).

## Usage

### Quickstart

If you have Docker and Nextflow installed, you can quickly test InterProScan and download the required data by running:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker,test \
  --datadir data \
  --interpro latest
```

Explanation of parameters:

* `-r 6.0.0-beta`: Specifies the version of InterProScan to run. We strongly recommend always specifying a version to ensure consistent and reproducible results.
* `-profile docker,test`:
  * `docker`: Executes tasks in Docker containers.
  * `test`: Uses a small test FASTA file included in the workflow.
* `--datadir data`: Sets the `data` directory as the location to store InterPro and member database files. The directory will be created automatically if it doesn't exist, and required files will be downloaded into it.
* `--interpro latest`: Uses the latest available InterPro data release.

> [!NOTE]
> While `--interpro latest` is the default, we strongly recommend pinning a specific version (e.g. `--interpro 106.0`) to ensure reproducibility.

After the run completes, the following files will be created in your working directory:

* `test.faa.gff3`: annotations in GFF3 format
* `test.faa.json`: Full annotations in JSON format
* `test.faa.tsv`: Tabular summary of matches (TSV format)
* `test.faa.xml`: Full annotations in XML format

The JSON and XML outputs are more comprehensive, the TSV is a concise summary, and the GFF3 is a standard format suitable for genome browsers and annotation pipelines.

### Scanning your own sequences

To annotate your own sequences FASTA file, omit the `test` profile and specify `--input`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa
```

For nucleotide sequences, add `--nucleic`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.fna \
  --nucleic
```

### Selecting specific analyses

To run only certain analyses, e.g. Pfam and MobiDB-lite:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa \
  --applications Pfam,MobiDB-lite
```

> [!TIP]
> Analysis names are case-insensitive and ignore hyphens and underscores: `MobiDB-lite`, `mobidblite`, and `MOBIDB_LITE` are all valid.

> [!NOTE]  
> Not sure which analyses to choose? See the [Available analyses](#available-analyses) table for descriptions of each tool and license requirements.

### Including GO terms and pathway annotations

To include Gene Ontology terms and pathway annotations:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa \
  --goterms \
  --pathways
```

### Running InterProScan on an HPC cluster with Slurm

To run InterProScan on your institute's Slurm cluster, use the `slurm` profile. This ensures that each task in the pipeline is submitted as a job to the Slurm scheduler.

Most HPC systems do not support Docker, but they often support Singularity or Apptainer for containerized execution. Include the appropriate profile (`singularity` or `apptainer`).

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile singularity,slurm \
  --datadir data \
  --input /path/to/sequences.faa
```

> [!IMPORTANT]
> The directory specified by `--datadir` must be accessible from all cluster nodes. This usually means it should be located on a shared network file system (e.g. NFS or Lustre).

## Available analyses

| Name             | Description                                                                             | Included |
|------------------|-----------------------------------------------------------------------------------------|----------|
| AntiFam          | Identifies sequences likely to be spurious or misannotated                              | ✅ Yes    |
| CATH-Gene3D      | Identifies structural domains from the CATH classification                              | ✅ Yes    |
| CATH-FunFam      | Groups protein domains into functional families based on CATH                           | ✅ Yes    |
| CDD              | Detects conserved domains using position-specific scoring matrices from NCBI            | ✅ Yes    |
| COILS            | Predicts coiled-coil regions based on sequence patterns                                 | ✅ Yes    |
| DeepTMHMM        | Predicts transmembrane helices                                                          | ❌ No     |
| HAMAP            | Identifies high-confidence protein families in microbial and organellar proteomes       | ✅Yes     |
| MobiDB-lite      | Predicts intrinsically disordered regions                                               | ✅ Yes    |
| NCBIFAM          | Matches proteins to curated HMMs from NCBI, including TIGRFAMs                          | ✅ Yes    |
| PANTHER          | Classifies proteins into families and subfamilies with curated GO terms                 | ✅ Yes    |
| Pfam             | Detects protein domains and families using HMMs built from multiple sequence alignments | ✅ Yes    |
| Phobius          | Predicts transmembrane topology and signal peptides                                     | ❌ No     |
| PIRSF            | Classifies proteins into evolutionary families based on full-length sequence similarity | ✅ Yes    |
| PIRSR            | Identifies conserved residues using manually curated site rules                         | ✅ Yes    |
| PRINTS           | Detects protein families using groups of conserved motifs                               | ✅ Yes    |
| PROSITE-patterns | Identifies protein features based on short sequence motifs                              | ✅ Yes    |
| PROSITE-profiles | Detects protein families and domains using position-specific scoring profiles           | ✅ Yes    |
| SFLD             | Classifies enzymes by relating sequence features to chemical function                   | ✅ Yes    |
| SMART            | Identifies signaling and extracellular domains                                          | ✅ Yes    |
| SUPERFAMILY      | Assigns structural domains using HMMs based on the SCOP superfamily classification.     | ✅ Yes    |
| SignalP-Euk      | Predicts signal peptides in eukaryotic proteins                                         | ❌ No     |
| SignalP-Prok     | Predicts signal peptides in prokaryotic proteins                                        | ❌ No     |

## Licensed analyses

DeepTMHMM, Phobius and SignalP contain licensed components and are disabled by default. 

> [!TIP]
> You do not need to install all three.
> Only download and configure the tool(s) you intend to use (e.g. just SignalP or Phobius).

To enable and execute any of these analyses:

1. Obtain a license for the tool.
2. Download and extract the archive.
3. Set the full path to the extracted directory in a Nextflow config file.

### Obtaining licensed components

#### DeepTMHMM 1.0

Request a standalone copy of DeepTMHMM 1.0 by sending an email to <licensing@biolib.com>. After receiving the package, extract it:

```sh
unzip -q DeepTMHMM-v1.0.zip
```

Then get the full path to the extracted directory:

```sh
echo "${PWD}/DeepTMHMM
```

#### Phobius 1.01

> [!IMPORTANT]  
> Phobius does not support certain non-standard or ambiguous residues. Any sequence containing pyrrolysine (one-letter code `O`), Asx (Asp/Asn ambiguity, `B`), Glx (Glu/Gln ambiguity, `Z`) or Xle (Leu/Ile ambiguity, `J`) will be skipped by Phobius but will continue to be processed normally by all other applications.

Download a copy of Phobius 1.01 [from Erik Sonnhammer's website](https://software.sbc.su.se/phobius.html), then extract:

```sh
tar -zxf phobius101_linux.tgz
```

And get the full path:

```sh
echo "${PWD}/phobius"
```

#### SignalP 6.0

SignalP 6.0 supports two modes: 

* Slow (full model)
* Fast (distilled model), recommended for most users

You need a license to download either model:

* full model: [slow_sequential](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=slow_sequential)
* distilled model: [fast](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast)

Extract the archive:

```sh
tar -zxf signalp-6.0h.fast.tar.gz
```

Then get the full path:

```sh
echo "${PWD}/signal6p_fast"
```

### Executing licensed analyses

You must define the tool path(s) in a Nextflow config file, such as `licensed.conf`.

#### Example: only Phobius

If you only want to run Phobius:

```groovy
params {
  appsConfig {
    phobius {
      dir = "/full/path/to/phobius"
    }
  }
}
```

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  -c licensed.conf \
  --input /path/to/sequences.faa \
  --applications phobius
```

#### Example: configuring multiple tools

To configure multiple licensed tools in one file:

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

And run with:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0-beta \
  -profile docker \
  -c licensed.conf \
  --input /path/to/sequences.faayour.fasta \
  --applications deeptmhmm,phobius,signalp_euk,signalp_prok \
  --no-matches-api
```

> [!WARNING]  
> As DeepTMHMM 1.0 and SignalP 6.0 predictions are not yet available in the [Matches API](https://www.ebi.ac.uk/interpro/matches/api/), the pre-calculated matches lookup must be disabled with `--no-matches-api`.

> [!NOTE]  
> Running both `signalp_euk` and `signalp_prok` will execute SignalP twice, once with eukaryotic post-processing and once without. Choose the mode best suited to your dataset.

## Documentation

Our full documentation is available on [ReadTheDocs](https://interproscan-docs.readthedocs.io/en/v6/).

## Support

For further assistance, please [create an issue](https://github.com/ebi-pf-team/interproscan6/issues) or [contact us](http://www.ebi.ac.uk/support/interproscan).

## Citation

If you use InterPro in your work, please cite the following publication:

> Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A. **InterPro: the protein sequence classification resource in 2025**. *Nucleic Acids Res*. 2025 Jan;53(D1):D444-D456. [doi: 10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).

