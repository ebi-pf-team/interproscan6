# InterProScan 6

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

> [!NOTE]  
> Phobius, SignalP and DeepTMHMM also require separate licenses and data downloads. See [Licensed analyses](#licensed-analyses) below.

## Usage

### Running the built‑in test

To test InterProScan, run the following command:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile test,docker \
  --datadir data \
  --interpro latest \
  --download
```

Explanation of parameters:

* `-profile test,docker`:
  * `test`: use an included example FASTA file
  * `docker`: execute tasks in Docker containers
* `--datadir data`: use `data` as the directory for storing all required databases; created automatically if needed
* `--interpro latest`: fetch the most recent InterPro release
* `--download`: download any missing metadata and database files
* `--max-workers`: Maximum number of workers available for the `InterProScan` when running locally

> [!IMPORTANT]
> *--max-workers* only applies when using the `local` profile (i.e. `-profile local`), it does **_not_** apply when running on a cluster.
> IPS6 will always use a minimum or 2 CPUs, with at least 1 dedicated to the main workflow and 1 to run
> processes (exception for PRINTS member, which require 2 CPUs to run processes).

After completion, you’ll find three output files in your working directory:

* `test.faa.json`: full annotations (JSON)
* `test.faa.tsv`: tabular summary (TSV)
* `test.faa.xml`: full annotations (XML)

The JSON and XML outputs are more comprehensive; the TSV is a concise summary.

### Scanning your own sequences

To annotate your own sequences FASTA file, omit the `test` profile and specify `--input`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile docker \
  --datadir data \
  --interpro latest \
  --input /path/to/sequences.fasta
```

For nucleotide sequences, add `--nucleic`:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile docker \
  --datadir data \
  --interpro latest \
  --input /path/to/sequences.fna \
  --nucleic
```

### Selecting specific analyses

To run only certain analyses (for example, Pfam and MobiDB‑lite):

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile test,docker \
  --datadir data \
  --interpro latest \
  --applications Pfam,MobiDB-lite
```

Analyses are case‑insensitive and ignore hyphens, so `MobiDB-lite` and `mobidblite` both work.

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

### Licensed analyses

DeepTMHMM, Phobius and SignalP contain licensed components and are disabled by default. To enable and execute these analyses, you first need to obtain their license, and data.

#### Obtaining licensed components

For each analyses, you need to request a license, then download and extract an archive that contains the data and machine learning model(s) required to run the analysis.

##### DeepTMHMM 1.0

Request a standalone copy of DeepTMHMM 1.0 by sending an email to <licensing@biolib.com>, then extract it when you receive it:

```sh
unzip -q DeepTMHMM-v1.0.zip
```

and make note of the full path to the package:

```sh
echo "${PWD}/DeepTMHMM
```

##### Phobius 1.01

Download a copy of Phobius 1.01 [from Erik Sonnhammer's website](https://software.sbc.su.se/phobius.html), then extract it:

```sh
tar -zxf phobius101_linux.tgz
```

and make note of the full path to the package:

```sh
echo "${PWD}/phobius"
```

##### SignalP 6.0

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

#### Executing licensed analyses

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
  -c licensed.conf \
  -profile docker \
  --input /path/to/sequences.fasta \
  --applications deeptmhmm,phobius,signalp_euk,signalp_prok \
  --offline
```

> [!WARNING]  
> DeepTMHMM 1.0 and SignalP 6.0 predictions are not yet available in the [Matches API](https://www.ebi.ac.uk/interpro/matches/api/). The pre-calculated matches lookup needs to be disabled with `--offline`.

> [!WARNING]  
> Phobius does not support certain non-standard or ambiguous residues. Any sequence containing pyrrolysine (one-letter code `O`), Asx (Asp/Asn ambiguity, `B`), Glx (Glu/Gln ambiguity, `Z`) or Xle (Leu/Ile ambiguity, `J`) will be skipped by Phobius but will continue to be processed normally by all other applications.

> [!NOTE]  
> Running both `signalp_euk` and `signalp_prok` will execute SignalP twice, once with eukaryotic post-processing and once without. Choose the mode best suited to your dataset.

## Documentation

Our full documentation is available on [ReadTheDocs](https://interproscan-docs.readthedocs.io/en/v6/).

# Citation

If you use InterPro in your work, please cite the following publication:

> Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A. **InterPro: the protein sequence classification resource in 2025**. *Nucleic Acids Res*. 2025 Jan;53(D1):D444-D456. [doi: 10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).

# Support

For further assistance, please [create an issue](https://github.com/ebi-pf-team/interproscan6/issues) or [contact us](http://www.ebi.ac.uk/support/interproscan).
