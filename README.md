# InterProScan6

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)
[![pytest](https://img.shields.io/badge/tested_with-pytest-337ab7.svg)]([https://github.com/askimed/nf-test](https://docs.pytest.org/en/stable/))
![Unit tests](https://github.com/ebi-pf-team/interproscan6/actions/workflows/unit-tests.yml/badge.svg)
[![codecov](https://codecov.io/gh/ebi-pf-team/interproscan6/graph/badge.svg?token=7MP9WCJHAQ)](https://codecov.io/gh/ebi-pf-team/interproscan6)  
[![Check Docker Files](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml/badge.svg)](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml)
[![Citation](https://github.com/ebi-pf-team/interproscan6/actions/workflows/citation.yml/badge.svg)](#citation)

> [!CAUTION]
> InterProScan6 is currently under active development and is not yet stable enough for a full release.

[InterPro](http://www.ebi.ac.uk/interpro/) is a database which integrates together predictive information about proteins’ function from a number of partner resources, giving an overview of the families that a protein belongs to and the domains and sites it contains.

Users who have novel nucleotide or protein sequences that they wish to functionally characterise can use the software package `InterProScan` to run the scanning algorithms from the InterPro database in an integrated way. Sequences are submitted in FASTA format. Matches are then calculated against all of the required member database’s signatures and the results are then output in a variety of formats.

# Documentation

Our full documentation is still under construction.

# Table of Contents
<!-- TOC -->
- [`InterProScan`](#interproscan6)
- [Documentation](#documentation)
- [Setup and Requirements](#set-up-and-requirements)
- [Requirements](#requirements)
- [Set up](#set-up)
  - [Quick installation](#quick)
  - [Install from source](#installing-from-source)
- [Using `InterProScan6`](#using-interproscan6)
  - [Quick start](#quick-start)
  - [Profiles](#profiles)
  - [Optional arguments](#optional-arguments)
  - [Applications](#applications)
  - [Using DNA sequences](#using-dna-sequences)
  - [Input sequences](#input-sequences)
  - [Outputs and results](#outputs)
- [Installing licensed applications (`MobiDB`, `Phobius`, `SignalP`, `TMHMM`)](#installing-licensed-applications-mobidb-phobius-signalp-tmhmm)
  - [DeepTMHMM](#deeptmhmm)
  - [MobiDB-Lite](#mobidb-lite)
  - [SignalP (version 6)](#signalp)
    - [Setting up SignalP](#set-up-1)
    - [Running SignalP](#running-interproscan6-with-signalp6)
    - [Changing mode](#changing-the-mode-of-signalp6-in-interproscan6)
    - [Converting from CPU to GPU](#converting-from-cpu-to-gpu-and-back-again)
- [Citing `InterProScan`](#citation)
- [Trouble shooting](#trouble-shooting)
<!-- /TOC -->

# Requirements

* `Java` (version >= 11)
* `Nextflow` (version >=23.04.02)
* A container run time
  * `InterProScan` includes built in support for:
    * `Docker` (version >= 24.0.5)
    * `SingularityCE` (version >= 4.2.0)
    * `Apptainer` (version >= 1.3.4)
  * Nextflow also supports using Charliecloud, Podman, Sarus, and Shifter. However you will need to build your own [Nextflow profiles](https://nextflow.io/docs/latest/config.html#config-profiles) to support these container runtimes

# Set up

The instructions below rely on an internet connection to pull the necessary images from Docker Hub.

## Quick

1. **Download InterPro data files**

```bash
# replace interpro-version with the appropriate version number
INTERPRO_VERSION="102.0"
curl "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/$INTERPRO_VERSION/interproscan-data-$INTERPRO_VERSION.tar.gz" \
    --output interproscan-data-<interpro-version>.tar.gz
tar -pxzf interproscan-data-<interpro-version>.tar.gz
mv interproscan-data-<interpro-version>/data .
rm interproscan-data-<interpro-version> -rf
rm interproscan-data-<interpro-version>.tar.gz
```

Running these commands within the `InterProScan` project directory should download and store all 
InterPro entry (XREF) and member database data in the `data` directory.

> [!IMPORTANT]
> By default `InterProScan` will look for a `data` directory in the `InterProScan` project dir. 
> If you store these data in an alternative directory use the `--datadir` flag to 
> provide the relative or absolute path to the directory when running `InterProScan`.

2. **Pull the `InterProScan6` base image from DockerHub using your container runtime of choice**

The base image includes all non-licensed dependencies including
`HMMER`, `BLAST`, `BioPython`, `easel`, etc.

Using `Docker`:
```bash
docker pull interproscan6:latest
```

Using `Singularity`:
```bash
singularity pull interproscan6.sif docker://interpro/interproscan6:latest
```

Using `Apptainer`:
```bash
apptainer pull interproscan6.sif docker://interpro/interproscan6:latest
```

> [!IMPORTANT]
> If you are using Singularity or Apptainer, `InterProScan6` expects to find the image files 
> within the current working directory. Otherwise the respective paths to the image files
> will need to be updated  in the respective `utilities/profiles/` config files.

3. **(Optional) Install licensed software**

By default `MobiDB`, `Phobius`, `SignalP`, and `TMHMM` member database analyses are deactivated in `InterProScan6` 
because they contain licensed components. In order to activate these analyses 
please see the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-tmhmm) documentation.

## Installing from source

1. **Download InterPro data file.**

```bash
# replace interpro-version with the appropriate version number
INTERPRO_VERSION="102.0"
curl "https://ftp.ebi.ac.uk/pub/databases/interpro/iprscan/6/$INTERPRO_VERSION/interproscan-data-$INTERPRO_VERSION.tar.gz" \
    --output interproscan-data-<interpro-version>.tar.gz
tar -pxzf interproscan-data-<interpro-version>.tar.gz
mv interproscan-data-<interpro-version>/data .
rm interproscan-data-<interpro-version> -rf
rm interproscan-data-<interpro-version>.tar.gz
```

2. **Build the Docker image.** (This includes the idrpred tool for MobiDB predictions). Run this command from the root of this repository.

```bash
docker build -t interproscan6 .
```

3. **(Optional) Install licensed software**

By default `MobiDB`, `Phobius`, `SignalP`, and `TMHMM` member database analyses are deactivated in `InterProScan6` 
because they contain licensed components. In order to activate these analyses 
please see the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-tmhmm) documentation.

4. **(Optional) Convert the Docker images to alternative container runtime images**

For example, to convert the ``interproscan6`` image to a Apptainer image, save the Docker image to a ``tar`` archive file then build an Apptainer image from the ``tar`` archive.

```bash
docker save interproscan6 > interproscan6.tar
apptainer build interproscan6.sif docker-daemon://interproscan6:latest
```

The same applies for Singularity:

```bash
docker save interproscan6 > interproscan6.tar
singularity build interproscan6.sif docker-daemon://interproscan6:latest
```

# Using `InterProScan6`

## Quick start

`InterProScan6` is configured via the command-line. The only mandatory arguments are the runtime profiles (`-profiles`) and input FASTA file (`--input`).

```bash
nextflow run interproscan.nf \
  -profile <container runtime, and executor> \
  --input <path to fasta file>
```

**`-profile` must** be included, and is used to define the executor and the container runtime used.

**`--input` must** be included, and is used to define the path to the input file containing the query sequences to be analysed in FASTA format.

> [!NOTE]  
> Note that `-profile` has a single dash because this is a Nextflow argument, but `--input` has two dashes because this is a `InterProScan6` argument.

## Profiles

A [Nextflow profile](https://nextflow.io/docs/latest/config.html#config-profiles) is a configuration file that defines runtime configuration attributes.

For `InterProScan6` to run, a profile for the container runtime and a profile for a executor must also be defined. `InterProScan6` includes the following profiles:

* Executor profiles:
  * local
  * slurm
  * lsf
* Container runtime profiles:
  * docker
  * singularity
  * apptainer

For example, to run `InterProScan` a cluster with the SLURM scheduler and Singularity:

```bash
nextflow run interproscan.nf \
  -profile slurm,singularity \
  --input <path to fasta file> 
```

Nextflow also supports using Charliecloud, Podman, Sarus, and Shifter. However you will need to build your own [Nextflow profiles](https://nextflow.io/docs/latest/config.html#config-profiles) to support these container runtimes.

## Optional arguments

**`--applications`** - Applications/member databases to run. By default `InterProScan` runs all member databases in the consortium ([except Mobidb-Lite due to licensing reasons](#mobidb)). Use the `--applications` to define a comma separate list of applications names (case insensitive).

**`--datadir`** - Path to the data directory. By default `InterProScan` looks for a `data` directory 
in the `InterProScan` project directory.

**`--disable_precalc`** - `InterProScan6` will check against a remote database of precalculated matches. The downstream analyses will then only be run against query sequences for whom no precalcualted match is available. You can disable this operation and run the analyses against all sequences in the input FASTA file by including the `--disable_precalc` flag in the `InterProScan6` command.

**`--formats`** - List output file formats as a comma separated list. Supported: JSON, TSV, XML. Default: TSV, JSON and XML

**`--goterms`** - Include Gene Ontology (GO) annotations in the final results.

**`--help`** - Display the help message.

**`--outdir`** - Output directory. By default the output files are written to the current working directory. Use `--outdir` to define the relative or abosolute path to the output directory. `InterProScan` will build all necessary directories.

**`--pathways`** - Include corresponding Pathway annotations in the final results. 

> [!WARNING]
> When using `--datadir` and `--outdir`, `InterProScan6` does not tolerate spaces in file paths.

For example, to run `InterProScan6` using only AntiFam and SFLD, without checking for pre-calculated matches in InterPro (using an example input file), with writing the results to the directory `results`, writing the results to a JSON and XML file and including GO term and Pathway annotation data in the final results, while using Docker on a local system:

```bash
nextflow run interproscan.nf \
  -profile docker,local \
  --input files_test/best_to_test.fasta \
  --applications signalp,antifam \
  --disable_precalc \
  --formats json,xml \
  --outdir results \
  --goterms \
  --pathways
```

> [!TIP]
> `InterproScan6` parameters are prefixed with a double dash, `--` (e.g. `--input`), Nextflow parameters are prefixed with a single dash, `-` (e.g. `-resume`)

## Applications

Below is a list of the applications (built in and those that require additional installation steps) that are available in `InterProScan6`:

* InterPro Consortium:
  * [Cath-Gene3D]( https://www.cathdb.info/) (use as 'Gene3D' to run Cath-Gene3D in `InterProScan`)
  * [CDD](https://www.ncbi.nlm.nih.gov/cdd)
  * [HAMAP](https://hamap.expasy.org/)
  * [MobiDB Lite](http://old.protein.bio.unipd.it/mobidblite/)
  * [NCBIfam](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/evidence/)
  * [PANTHER](http://www.pantherdb.org/)
  * [Pfam](https://pfam.xfam.org/)
  * [PIRSF](https://proteininformationresource.org/pirsf/)
  * [PRINTS](https://interpro-documentation.readthedocs.io/en/latest/prints.html)
  * [PROSITE profiles](https://prosite.expasy.org/)
  * [SFLD](http://sfld.rbvi.ucsf.edu/archive/django/index.html)
  * [SMART](http://smart.embl-heidelberg.de/)
  * [SUPERFAMILY](https://supfam.mrc-lmb.cam.ac.uk/)
* Applications not in the consortium:
  * Antifam
  * [COILS](http://www.ch.embnet.org/software/COILS_form.html)
  * [DeepTMHMM](https://www.biorxiv.org/content/10.1101/2022.04.08.487609v1) (*not yet available*)
  * [FunFam](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-019-2988-x)
  * [Phobius](http://phobius.sbc.su.se/)
  * [PIRSR](https://www.uniprot.org/help/pir_rules)
  * [PROSITE patterns](http://prosite.expasy.org/)
  * [SignalP](http://www.cbs.dtu.dk/services/SignalP/) (sets the organism argument for SignalP6 to `other`)
  * SignalP_EUK (sets the organism argument for SignalP6 to `eukaryote`)

> [!NOTE]
> Quoting the SignalP documentation:  
> Specifying the eukarya method of `SignalP6` (`SignalP_EUK`) triggers post-processing of the SP predictions by `SignalP6` to prevent spurious results (it will only predicts type Sec/SPI).

## Using DNA sequences

`InterProScan6` takes advantage of the Open Reading Frame (ORF) prediction tool `esl-translate` within the [`easel` tool suite](https://github.com/EddyRivasLab/easel) in order to analyse nucleic acid sequences.

The `easel` application itself and all of its dependencies are integrated in InterProScan via a docker image (see [set-up](#set-up)).

**To run anlyses with nucleic acid sequences, run `InterProScan6` with the `--nucleic` flag**

    nextflow run interproscan.nf \
        --input <path to fasta file> \
        -profile <executor,container runtime> \
        --nucleic

By default `InterProScan6` will assume the input FASTA file contains protein sequences. The `--nucleic` flag instructs `InterProScan6` to retrieve all possible ORFs using the `easel` tool suite.

> [!WARNING]  
> The input FASTA file must contain sequences of the same type, i.e. _all_ protein sequences or _all_ nucleic acid sequences.

You can configure the prediction of ORFs by updating the relevant `translate` parameters in `nextflow.config`:

```groovy
    translate { 
        strand = 'both'  
        methionine = false  
        min_len = 20
        genetic_code = 1
    }
```

* `strand` - DNA strand(s) to be translated
    - `'both'`
    - `'plus'`
    - `'minus'`
* `methionine` - predicted ORFs start with M (methionine)
    - `false` - use initation codon
    - `true` - all ORFs start with M
* `min_len` - minimum length of predicted ORFs [any interger]
* `genetic_code` - ID of the genetic code to use

<table>
  <thead>
    <tr>
      <th>ID</th>
      <th>Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>1</td>
      <td>Standard</td>
    </tr>
    <tr>
      <td>2</td>
      <td>Vertebrate mitochondrial</td>
    </tr>
    <tr>
      <td>3</td>
      <td>Yeast mitochondrial</td>
    </tr>
    <tr>
      <td>4</td>
      <td>Mold, protozoan, coelenterate mitochondrial; Mycoplasma/Spiroplasma</td>
    </tr>
    <tr>
      <td>5</td>
      <td>Invertebrate mitochondrial</td>
    </tr>
    <tr>
      <td>6</td>
      <td>Ciliate, dasycladacean, Hexamita nuclear</td>
    </tr>
    <tr>
      <td>9</td>
      <td>Echinoderm and flatworm mitochondrial</td>
    </tr>
    <tr>
      <td>10</td>
      <td>Euplotid nuclear</td>
    </tr>
    <tr>
      <td>11</td>
      <td>Bacterial, archaeal; and plant plastid</td>
    </tr>
    <tr>
      <td>12</td>
      <td>Alternative yeast</td>
    </tr>
    <tr>
      <td>13</td>
      <td>Ascidian mitochondrial</td>
    </tr>
    <tr>
      <td>14</td>
      <td>Alternative flatworm mitochondrial</td>
    </tr>
    <tr>
      <td>16</td>
      <td>Chlorophycean mitochondrial</td>
    </tr>
    <tr>
      <td>21</td>
      <td>Trematode mitochondrial</td>
    </tr>
    <tr>
      <td>22</td>
      <td>Scenedesmus obliquus mitochondrial</td>
    </tr>
    <tr>
      <td>23</td>
      <td>Thraustochytrium mitochondrial</td>
    </tr>
    <tr>
      <td>24</td>
      <td>Pterobranchia mitochondrial</td>
    </tr>
    <tr>
      <td>25</td>
      <td>Candidate Division SR1 and Gracilibacteria</td>
    </tr>
  </tbody>
</table>

## Input sequences

* The input file must be in FASTA format
* The sequences in the input file must be _all_ protein or _all_ nucleic acid sequences
* The table below lists the illegal characters for each member db with illegal characters:

<table>
  <thead>
  <tr>
    <th>Member DB</th>
    <th>Illegal characters</th>
  </tr>
  </thead>
  <tbody>
  <tr>
    <td>AntiFam</td>
    <td>-</td>
  </tr>
  <tr>
    <td>FunFam</td>
    <td>- _ .</td>
  </tr>
  <tr>
    <td>Gene3D</td>
    <td>- _ .</td>
  </tr>
  <tr>
    <td>HAMAP</td>
    <td>- _ .</td>
  </tr>
  <tr>
    <td>NCBIFAM</td>
    <td>-</td>
  </tr>
  <tr>
    <td>Panther</td>
    <td>-</td>
  </tr>
  <tr>
    <td>Pfam</td>
    <td>-</td>
  </tr>
  <tr>
    <td>Phobius</td>
    <td>- _ . * o x u z j</td>
  </tr>
  <tr>
    <td>PIRSR</td>
    <td>-</td>
  </tr>
  <tr>
    <td>PIRSF</td>
    <td>-</td>
  </tr>
  <tr>
    <td>PROSITE Profiles</td>
    <td>- _ .</td>
  </tr>
  <tr>
    <td>SFLD</td>
    <td>- _ .</td>
  </tr>
  <tr>
    <td>SUPERFAMILY</td>
    <td>-</td>
  </tr>
  </tbody>
</table>

## Outputs

The output can be retrieved in TSV, XML or JSON format.

### TSV

The TSV format provides a summary of **only** the matches found. Each row represents a unique match between a signature and a query sequence.

1.  Protein accession
2.  Sequence MD5 digest
3.  Sequence length
4.  Source Member database for the signature
5.  Signature accession
6.  Signature description
7.  Match start location in the query sequence
8.  Match end location in the query sequence
9.  Score of the match
    * E-value for AntiFam, Cath-Gene3D, FunFam, NCBIFam, PANTHER, Pfam, PIRSF, PRINTS, SFLD, SMART, SUPERFAMILY, CDD
    * Score for HAMAP, PROSITE ProSiteProfiles
    * '-' for Coils, MobiDB-lite, Phobius, PROSITE Patterns, SignalP, TMHMM
10. Status of the match (T: true)
11. Date of the run (format `DD-MM-YYYY`)
12. Accession of the associataed InterPro entry
13. Description from the Associated InterPro entry
14. Pipe-separated list of GO annotations with their source(s). Only displayed if the `--goterms` option is switched on
15. Pipe-separated list of pathways annotations. Only displayed if the `--pathways` option is switched on

### JSON

The JSON file contains more data than the TSV, and lists information for _all_ proteins in the input FASTA file, not only those for whom matches were found.

For each input/query sequence:

1. `sequence`: The submitted protein or nucleotie sequence
2. `md5`: MD5 hash of the submitted sequence
3. `matches`: List of matches (JSON objects). For each match:
    * `evalue`: Overall, full sequence evalue
    * `score`: Overall, full sequence bit-score
    * `model-ac`: Accession of the member database model
    * `signature`: A JSON object summarising the InterPro signature
        * `accession`: Signature accession
        * `name`: Name from the InterPro entry
        * `description`: Description from the InterPro entry
        * `entry`: The accession of the InterPro entry that the signature is associated with
            * entry will be null if a singature is not associated with an InterPro entry
            * `accession`: The InterPro entry accession
            * `name`: The InterPro entry name
            * `description`: The InterPro entry description
            * `type`: The type of InterPro entry (e.g. family, domain, etc.)
            * `goXRefs`: Geneontology (GO) terms associated with the InterPro entry - only retrieved if the `--goterms` flag is used
            * `pathwayXRefs`: Pathway information associated with the InterPro entry - only retrieved if the `--pathways` flag is used
        * `signatureLibraryRelease`: JSON object containing:
            * `library`: Application/member database name
            * `version`: Release version number
    * `locations` : List of locations where the signature matched the protein sequence. Specifically, this is a list of JSON objects, one JSON object per location where a match between the protein sequence and signature was found. For each location:
        * `start` Start point of the alignment location with respect to the query sequence
        * `end` End point of the alignment location with respect to the query sequence
        * `hmmStart` Start point of the local alignment with respect to the HMM profile
        * `hmmEnd` End point of the local alignment with respect to the HMM profile
        * `evalue`: Independent E-value
        * `score`: Bit score
        * `envelopesStart`: Start of the envelop
        * `envelopeEnd`: End of the envelop
        * `location-fragments`: List of JSON objects, one JSON object per fragment:
            * `start`: Start location of the fragment in the query sequence
            * `end`: End location of the fragment in the query sequence
            * `dc-status`: Continuous/discontinuous status.
        * `sites`: List of JSON objects, one JSON object per site (a domain signature in some member databases can have multiple sites). Per site:
            * `description`: Site description (from InterPro)
            * `numLocations`: The number of locations
            * `siteLocations`: List, one JSON object per location:
                * `start`: Start location of the site in the query sequence
                * `end`: End location of the site in the query sequence
                * `residue`: The amino acid residue of the site
    * `xref`: The protein sequence ID and description listed in the input FASTA file

### XML

The richest form of the data is the XML representtaion, and includes data for all sequences 
listed in the input FASTA File.

The XML Schema Definition (XSD) is available [here](http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/schemas/). `InterProScan6` uses the latest XSD.

For each query sequence:

1. `sequence`: The submitted protein or nucleotie sequence
2. `xref`: The sequence ID and name/description from the input FASTA file
3. `md5`: MD5 hash of the submitted sequence
4. `matches`: List of matches:
    * `<>-match`: The type/source of the signature match:
      * `hmmer3-match`: AntiFam, NCBIFam, FunFam, Gene3D, HAMAP, Panther, PFAM, PIRSF, PIRSR, SFLD, SMART, SUPERFAMILY
      * `<member-name>-match`: Match from a member database that does not use HMMER, e.g. CDD
    * The information for both these keys is very similar and is summarised here:
        * `signature`: Represents the member database signature. Includes accession, name and description
            * `entry`: Associated InterPro entry. Includes entry accession, description, name and type (e.g. family, domain, etc.), as well as any associated pathway information (if the `--pathway` flag is used) and Geneontology (GO) terms (of the `--goterms` flag is used)
            * `library release`: Release version of the member datbase. Includes name and version/release number
        * `models`: Information about the model, including the name, and accession
        * `locations`: Represents domain hits in the query sequence. Includes:
            * E-value
            * Score: The bitscore or other member database relevant score
            * The envelop start and end: Start and end point of the envelop
            * Hmm-start and hmm-end: Start and end point of the local alignemnt with respect to the HMM profile
            * Hmm-length: Length of the alignemnt along the query sequence
            * Hmm-bounds: Description of the HMMER Hmm bound pattern
            * start and end: Start and end point of the alignment location with respect to the query sequence
            * alignemnt: The query sequence alignment to the model
            * cigar-alignemnt: The [cigar alignment](https://replicongenetics.com/cigar-strings-explained/)
            * `site-loctaions`: information about sites (for those member databases that contain site data):
                * Each site is represented by a `site-location`, which has a start, stop and residue.

### Continuous and discontinuous status

The `dc-status` refers to continuous nature of a domain 
hit in some member databases. If a domain is not  continuous, i.e. 
is broken up into fragments, each of the fragments are represented under the
`location-fragments` key and are labelled as:
* "C_TERMINAL_DISC" - the most c-terminal fragment
* "N_TERMINAL_DISC" - the most n-terminal fragment
* "NC_TERMINAL_DISC" - all other fragments

> [!NOTE]  
> Not all member database analyses can detect discontinious domains. At the present only Gene3D and 
> FunFam are able to detect discontinious domains.

### The envelop

The envelope represents the region of a protein sequence where the domain may be located. Often it is wider than what HMMER chooses as a reasonably confident alignment.

**Panther exception:** The output from HMMER3 against the HMM models of Panther is post-processed to select only the best homologous family. Therefore, there is a maximum of one domain hit for each Panther signature in a protein sequence. Owing to this the E-value and Score and listed under the `signature` key, not the `locations` key.

# Installing licensed applications (`MobiDB`, `Phobius`, `SignalP`, `TMHMM`)

By default `MobiDB`, `Phobius`, `SignalP`, and `DeepTMHMM` member database analyses are deactivated in `InterProScan6` because they contain licensed components. In order to activate these analyses please obtain the relevant licenses and files from the provider (ensuring the software version numbers are the same as those supported by your current `InterProScan6` installation).

Files can be placed in any location.

> [!NOTE]  
> The instructions below presume `InterProScan6` is being run with the Docker. If an alterantive 
> container runtime is used the methods below will need to be adapated accordingly.
> As above, if using Singularity and Apptainer, the images should be kept in the root of the 
> `InterProScan6` directory. Otherwise please update the container paths in the respective `utilities/profiles/` config files.

## DeepTMHMM

Coming soon...

## MobiDB-Lite

Some of the compoments within `MobiDBLite` are GPL-licensed, meaning all software and data, and thus 
work that uses this software, also needs to be GPL-licensed. This may not be ideal or suitable
for all users. Therefore, we provide a version of the `MobiDBLite` analytical software that 
is not GPL-licensed, called [`idrpred`](https://github.com/matthiasblum/idrpred).

To setup `MobiDB`/`idrpred` for `InterProScan6` pull the `idrpred` Docker image from Docker hub using your container runtime of choice.

Using docker:
```bash
docker pull idrpred:latest
```

Using `Singularity`:
```bash
singularity pull idrpred.sif docker://matblum/idrpred/idrpred:latest
```

Using `Apptainer`:
```bash
apptainer pull idrpred.sif docker://matblum/idrpred/idrpred:latest
```

## `Phobius`

1. Download Phobius from the [Phobius server](https://software.sbc.su.se/phobius.html)

2. Unpack the `tar` file
```bash
tar -xzf phobius101_linux.tgz -C <PHOBIUS-DIR>
```

3. Copy the docker file available in the `./docker_files/phobius/` directory to your local `Phobius` directory
```bash
# with the terminal pointed at the root of this repo
cp docker_files/phobius/Dockerfile <PHOBIUS-DIR>/Dockerfile
```

4. Build a docker image -
```bash
# with the terminal pointed at your local phobius dir
docker image build -t phobius .
```

5. Check the `subworkflows/sequence_analysis/members.config` file to make sure the `Phobius` version is correct.
```groovy
    phobius {
            release = "1.01" <---- update if necessary
            runner = "phobius"
        }
```

6. (Optional) Convert the Docker image to an image of your container runtime.

For example, to build a singularity image:
```bash
docker save phobius > phobius.tar
singularity build phobius.sif docker-archive://phobius.tar
```

## `SignalP`

### Set up

1. Obtain a license and download `SignalP6` (`SignalP` version 6) from the [SignalP6 server](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (under 'Downloads').
    * Either fast or slow models can be implemented
    * To change the implemented mode please see the [Changing mode](#changing-mode) documentation

2. Unpackage the `SignalP6` `tar` file

```bash
tar -xzf signalp-6.0h.fast.tar.gz -C <SIGNALP-DIR>
```

3. Copy the docker file available in the `./docker_files/signalp/` directory to your local `SignalP6` directory

```bash
# with the terminal point at the root of this repo
cp docker_files/signalp/Dockerfile <SIGNALP-DIR>/Dockerfile
```

4. Build a docker image - _the Nextflow pipeline needs all third party tools to be stored within linux containers_. 

```bash
# with the terminal pointed at your local signalp dir
docker build -t signalp6 .
```

5. Check the version number in `subworkflows/sequence_analysis/members.config` is correct:

```groovy
signalp {
    release = "6.0h"  <--- make sure the release is correct
    runner = "signalp"
    ...
}
...
signalp_euk {
    release = "6.0h"  <--- make sure the release is correct
    runner = "signalp_euk"
    ...
}
```

6. (Optional) Convert the Docker image to an image of your container runtime.

For example, to build a singularity image:
```bash
docker save signalp6 > signalp6.tar
singularity build signalp6.sif docker-archive://signalp6.tar
```

### Running `InterProScan6` with `SignalP6`

Include `signalp` or `signalp_euk` in the list of applications defined using the `--applications` flag.

```bash
nextflow run interproscan.nf \
    --input utilities/test_files/best_to_test.fasta \
    --applications signalp \
    -profile local,docker
```

```bash
nextflow run interproscan.nf \
    --input utilities/test_files/best_to_test.fasta \
    --applications signalp_euk \
    -profile local,docker
```

> [!NOTE]
> Quoting the SignalP documentation:  
> Specifying the eukarya method of `SignalP6` (`SignalP_EUK`) triggers post-processing of the SP predictions by `SignalP6` to prevent spurious results (it will only predicts type Sec/SPI).

### Changing the mode of `Signalp6` in `InterProScan6`

`SignalP6` supports 3 modes: `fast`, `slow` and `slow-sequential`. The mode can be set using the `--signalp_mode` flag. The default mode is `fast`.

> [!WARNING]
> The slow mode can take 6x longer to compute. Use when accurate region borders are needed.

You may need to install the other models mannually, please see the [SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#installing-additional-modes).

For example, to run `InterProScan` with the input file `best_to_test.fasta`, using SignalP with only eukaryotic models in slow mode, and with retrieving precalculated matches disabled on a local machine using docker:

```bash
nextflow run interproscan.nf \
  --input utilities/test_files/best_to_test.fasta \
  --applications signalp_euk \
  --disable_precalc \
  --signalp_mode slow \
  -profile docker,local
```

> [!NOTE]  
> `InterProScan6` only supports implementing one `SignalP` mode at a time.

### Converting from CPU to GPU, and back again

By default, `SignalP` runs on your CPU. If you have a GPU available, you can convert the `SignalP` models so that your installation can use GPU-acceleration.

You will need to install `SignalP` in order to convert to GPU models.

1. Convert the ``SignalP`` installation to GPU by following the [SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#converting-to-gpu)

2. Build a docker image for `SignalP` with GPU

```bash
# with the terminal pointed at your local signalp dir
docker image build -t signalp6_gpu .
```

3. (Optional) Convert the image to your container runtime of choice

For example, to build a singularity image:
```bash
docker save signalp6_gpuu > signalp6_gpu.tar
singularity build signalp6_gpu.sif docker-archive://signalp6_gpu.tar
```

To run `SignalP` with GPU acceleration with `InterProScan6` use the flag `--signalp_gpu`.

For example, to run ``InterProScan`` with only ``SignalP`` enabled, using GPU acceleration on a SLURM cluster with Singularity support:

```bash
nextflow run interproscan.nf \\
  --input <fasta file> \\
  --applications signalp \\
  --signalp_gpu \\
  -profile singularity,slurm
```

# Benchmarking and trouble shooting the performance

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating a HTML report. However, these reports are limited to presenting the resource usage from only a single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to enable assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size, the number of CPUs or amount of memory allocated. 

## Using `Nextflow`

### Report

Run `IPS6` with the `-with-report [file name]` Nextflow parameter to generate a summary html report for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#execution-report).

### Timeline

Run `IPS6` with the `-with-timeline [file name]` Nextflow parameter to generate a visual representation of the operational timeline for the current run. You can find more in the [Nextflow docs](https://www.nextflow.io/docs/latest/tracing.html#timeline-report).

## Using IPS6 benchmarking scripts

1. **Setup**

Install all necessary third party packages listed in `benchmarking/requirements.txt`.

2. **Build a trace file as part of your InterProScan runs**

Include the following code in the `nextflow.config` file and rename the resulting trace file
`ips6.trace.txt` after each run.

```groovy
trace {
    enabled = true
    fields = 'task_id,process,realtime,duration,status,cpus,%cpu,time,memory,%mem,rss,peak_rss,submit,start,complete,queue'
    file = 'ips6.trace.txt'
    overwrite = true   
}
```

You can find more information the columns that can be included in the trace in the [Nextflow documentation](https://www.nextflow.io/docs/latest/tracing.html#trace-report). 

At a minimum, the `IPS6` benchmarking scripts require the following columns:
* process
* realtime
* rss
* peak_rss

3. **Update the benchmarking config file**

Update the `tracefile.json` file (or create a new JSON file) keyed by the name of each group of analyses,
e.g. the number of CPU assigned to the group or the batch size used, and valued by a list of string 
representations of paths to `IPS6` trace files.

For example, if you wanted to assess the impact of altering the batch size on performance:

```json
{
    "500": [
        "benchmarking/24.08.27.batch.500.report.3.tsv",
        "benchmarking/24.08.26.batch.500.report.2.tsv",
        "benchmarking/24.08.26.batch.500.report.1.tsv"
    ],
    "1000": [
        "benchmarking/24.08.27.batch.1000.report.3.tsv",
        "benchmarking/24.08.26.batch.1000.report.2.tsv",
        "benchmarking/24.08.26.batch.1000.report.1.tsv"
    ],
    "5000": [
        "benchmarking/24.08.27.batch.5000.report.3.tsv",
        "benchmarking/24.08.26.batch.5000.report.2.tsv",
        "benchmarking/24.08.26.batch.5000.report.1.tsv"
    ]
}
```

4. **Run the benchmarking**

The onl required argument for running the benchmarking is the path to the JSON file listing the 
paths to the trace files.

```bash
# running from the root of the IPS6 project dir
# and using the benchmarking/tracefiles.json file
python3 benchmarking/benchmark_ips6.py benchmarking/tracefiles.json
```

You can print a help message to prin out the argument options:
```bash
python3 benchmarking/benchmark_ips6.py --help
```

By default, the benchmarking will label the groupings as 'Groups' on the resulting plot axes and 
legends. You can name the groupings using the `--group_name` flag and providing the name you 
wish to be assigned to the axes and legends, e.g. `--group_name "Batch Sizes"`, or `--group_name "Number of CPU"`.

By default, the resulting figures are only written out in `PDF` format. Use the `--format` flag to 
list the desired file outputs. Accepted outputs: png, pdf, and svg. For example to generate svg and 
png files use `--format png,svg`.

By default the trace file writes the in human readable format, but can be configured to write the raw
values. If this is the case, include the `--raw` flag in the `benchmark_ips6.py` command.

By default, the output figures will be written to the current working directory. To write the files 
to a desired output directory use the `--outdir` flag and provide the path for the output dir. The 
scripts will build all necessary parent directories for the output dir.

If you wish to perform further analyses on the data, use the `--save_data` flag to configure 
`benchmark_ips6.py` to write out the dataframe it generates to a CSV file in the output dir.

## Output:

Each run of `benchmark_ips6.py` will produce the following figures (note, references to 'group' refers to the keys in the input JSON file, each key represents a different 'group'):

1. `total_runtime.*` - Shows the total run time of IPS6 per group in the input JSON file
2. `process_runtime.*` - Shows the total run time per process in IPS6
3. `process_runtime_piechart.*` - Shows the percentage of the total runtime contributed by each process
4. `pie_chart_values.csv` - Contains the data used to build the `process_runtime_piechart.*` figure. If many processes are included the legends in the pie chart can often overlap. Use this CSV file to plot the pie chart (or alternative chart).
5. `overall_memory_usage.*` - Plots the overall memory usage per group in the input JSON file
6. `overall_max_memory_usage.*` - Plots the overall maximum memory used per group in the input JSON file
7. `memory_per_process.*` - Plots the memory usage per process (and per group if multiple groups are defined in the input JSON file)
8. `max_memory_per_process.*` - Plots the maximum memory usage per process (and per group if multiple groups are defined in the input JSON file)

# Citation

If you use `InterProScan` or `InterPro` in your work, please cite the following publications:

**`InterProScan`:**

> Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H, Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A, Scheremetjew M, Yong SY, Lopez R, Hunter S. InterProScan 5: genome-scale protein function classification. Bioinformatics. 2014 May 1;30(9):1236-40. doi: 10.1093/bioinformatics/btu031. Epub 2014 Jan 21. PMID: 24451626; PMCID: PMC3998142.

**InterPro:**

> Paysan-Lafosse T, Blum M, Chuguransky S, Grego T, Pinto BL, Salazar GA, Bileschi ML, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunić I, Marchler-Bauer A, Mi H, Natale DA, Orengo CA, Pandurangan AP, Rivoire C, Sigrist CJA, Sillitoe I, Thanki N, Thomas PD, Tosatto SCE, Wu CH, Bateman A. InterPro in 2022. Nucleic Acids Res. 2023 Jan 6;51(D1):D418-D427. doi: 10.1093/nar/gkac993. PMID: 36350672; PMCID: PMC9825450.

You can find more information about the InterPro consortium (including the member databases) on the [InterPro website](https://www.ebi.ac.uk/interpro/about/consortium/).

# Trouble shooting

## Permission denied

On some systems, Nextflow requires root privileges to be able to create the output directories. 

If you receive an error message such as:

```bash
ERROR ~ Error executing process > 'PARSE_SEQUENCE (1)'

Caused by:
  Unable to create directory=.../InterProScan6/work/56/c74d6bbc6637c201a68afd84b75d27 -- check file system permissions
```

Try running Nextflow with root privileges:

```bash
sudo nextflow run interproscan.nf --input <path to fasta file> 
```

Also try providing root privileges to docker within Nextflow, by changing the the `runOptions` key in `nextflow.config`:

```
 docker {
     enabled = true
     mountFlags = 'Z'
+    runOptions = '--user root'  <<-- ADD
 }
```

## File not found

If you recieve a file not found error:

```bash
FileNotFoundError: [Errno 2] No such file or directory
```

This may be due to the docker image not being built correctly. This can happen due to file permissions from preventing docker from reading the `InterProScan6` files.

Try running docker with root privileges:

```bash
sudo docker build -t interproscan6 .
```

Check the docker installtion is configured correctly, with all necessary privileges. [StackOverflow](https://stackoverflow.com/questions/48957195/how-to-fix-docker-got-permission-denied-issue)

For example, check root privileges have been provided to the docker socket

```bash
sudo chmod 666 /var/run/docker.sock
```

## Segmentation fault

If you a recieve an error such as the following:

```bash
Command error:
  .command.sh: line 2:     7 Segmentation fault      (core dumped) /opt/hmmer3/bin/hmmsearch --cut_ga --cpu 1 -o 7.0._.antifam._.out AntiFam.hmm mini_test.1.fasta
```

This is generally due to HMMER being unable to find a necessary data file.
Make sure the data directory is correctly structured and populated and `InterProScan` is 
pointed to the correct data directory using the `--data` flag if not using the default data
directory location in the project dir.

## Segmentation fault

If you a recieve an error such as the following:
```bash
Command error:
  .command.sh: line 2:     7 Segmentation fault      (core dumped) /opt/hmmer3/bin/hmmsearch --cut_ga --cpu 1 -o 7.0._.antifam._.out AntiFam.hmm mini_test.1.fasta
```
This is generally due to HMMER being unable to find a necessary data file.
Make sure the data directory is correctly structured and populated and `InterProScan` is 
pointed to the correct data directory using the `--data` flag if not using the default data
directory location in the project dir.

## Cannot access or failed to open output files for writing

For example:

```bash
Command error:
  
  Error: Failed to open output file hmmer_AntiFam.hmm.out for writing
```

This is most likely a file permission error.

A potential fix is to provide root privilges to the docker contains run by Nextflow in `nextflow.config`:

```groovy
process.container = 'interproscan6'
docker {
    enabled = true
    runOptions = '--user root'
}
```
