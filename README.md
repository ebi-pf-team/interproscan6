# InterProScan 6

[![Nextflow](https://img.shields.io/badge/%E2%89%A525.04.6-0dc09d?style=flat&logo=nextflow&label=nextflow&labelColor=f5fafe)](https://www.nextflow.io/)
[![run with docker](https://img.shields.io/badge/docker-1D63ED?logo=docker&logoColor=1D63ED&label=run%20with&labelColor=f5fafe)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/singularity-1E4383?label=run%20with&labelColor=f5fafe)](https://sylabs.io/docs/)

[InterPro](http://www.ebi.ac.uk/interpro/) is a database that brings together predictive information on protein function from multiple partner resources. It provides an integrated view of the families, domains and functional sites to which a given protein belongs.

**InterProScan** is the command‑line tool that allows you to scan protein or nucleotide sequences against the InterPro member‑database signatures in a single workflow. Researchers with novel sequences can use InterProScan to annotate their data with family classifications, domain architectures and site predictions.

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
> Some InterProScan analyses rely on third-party applications that require separate licenses and downloads. This does not mean you need a license to run InterProScan itself or to analyse your own sequences with the bundled analyses.

## Usage

### Quickstart

If you have Docker and Nextflow installed, you can quickly test InterProScan and download the required data by running:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker,test \
  --datadir data \
  --interpro latest
```

Explanation of parameters:

* `-r 6.0.0`: Specifies the version of InterProScan to run. We strongly recommend always specifying a version to ensure consistent and reproducible results.
* `-profile docker,test`:
  * `docker`: Executes tasks in Docker containers.
  * `test`: Uses a small test FASTA file included in the workflow.
* `--datadir data`: Sets the `data` directory as the location to store InterPro and member database files. The directory will be created automatically if it doesn't exist, and required files will be downloaded into it.
* `--interpro latest`: Uses the latest available InterPro data release.

> [!NOTE]
> While `--interpro latest` is the default, we strongly recommend pinning a specific version (e.g. `--interpro 107.0`) to ensure reproducibility.

After the run completes, the following files will be created in your working directory:

* `test.faa.gff3`: Annotations in GFF3 format
* `test.faa.json`: Full annotations in JSON format
* `test.faa.jsonl`: Full annotations in JSON Lines format (one line for each input sequence)
* `test.faa.tsv`: Tabular summary of matches (TSV format)
* `test.faa.xml`: Full annotations in XML format

The JSON, JSON Lines, and XML outputs are more comprehensive, the TSV is a concise summary, and the GFF3 is a standard format suitable for genome browsers and annotation pipelines.

### Scanning your own sequences

To annotate your own sequences FASTA file, omit the `test` profile and specify `--input`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa
```

For nucleotide sequences, add `--nucleic`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.fna \
  --nucleic
```

## Documentation

Our full documentation is available on [ReadTheDocs](https://interproscan.readthedocs.io/).

## Support

For further assistance, please [create an issue](https://github.com/ebi-pf-team/interproscan6/issues) or [contact us](http://www.ebi.ac.uk/support/interproscan).

## Citation

If you use InterProScan in your work, please cite the following publication:

> Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A. **InterPro: the protein sequence classification resource in 2025**. *Nucleic Acids Res*. 2025 Jan;53(D1):D444-D456. [doi: 10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).
