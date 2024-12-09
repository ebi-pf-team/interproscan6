# InterProScan6

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)
![Unit tests](https://github.com/ebi-pf-team/interproscan6/actions/workflows/unit-tests.yml/badge.svg)
[![Check Docker Files](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml/badge.svg)](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml)
[![Citation](https://github.com/ebi-pf-team/interproscan6/actions/workflows/citation.yml/badge.svg)](#citation)

> [!CAUTION]
> InterProScan6 is currently under active development and is not yet stable enough for a full release.

[InterPro](http://www.ebi.ac.uk/interpro/) is a database which integrates together predictive information about proteins’ functions from a number of partner resources, giving an overview of the families that a protein belongs to as well as the domains and sites it contains.

Users who have novel nucleotide or protein sequences that they wish to functionally characterise can use the software package `InterProScan` to run the scanning algorithms from the InterPro database in an integrated way. Sequences are submitted in FASTA format. Matches are then calculated against all the required member databases signatures and the results are then output in a variety of formats.

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
- [Installing licensed applications (`Phobius`, `SignalP`, `TMHMM`)](#installing-licensed-applications-phobius-signalp-tmhmm)
  - [DeepTMHMM](#deeptmhmm)
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
* `Nextflow` (version >=23.04.01)
* A container run time
  * `InterProScan` includes built in support for:
    * `Docker` (version >= 24.0.5)
    * `SingularityCE` (version >= 4.2.0)
    * `Apptainer` (version >= 1.3.4)
  * Nextflow also supports using Charliecloud, Podman, Sarus, and Shifter. However you will need to build your own [Nextflow profiles](https://nextflow.io/docs/latest/config.html#config-profiles) to support these container runtimes

> [!WARNING]
> Support for Java versions prior to 17 were dropped in [Nextflow version 24.11.0-edge](https://www.nextflow.io/docs/latest/install.html).

# Set up

The instructions below rely on an internet connection to pull the necessary images from Docker Hub.

> [!IMPORTANT]
> By default `InterProScan` will look for a `data` directory in the `InterProScan` project dir.
> If you store these data in an alternative directory use the `--datadir` flag to
> provide the relative or absolute path to the directory when running `InterProScan`.

## Quick

Instead of setting up an installation, Nextflow can pull down the latest `InterProScan` image from DockerHub if one is not already available on the system being used to run `InterProScan`, otherwise Nextflow will use the available `InterProScan` image.

1. **Download InterPro data files**

```bash
curl -OJ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/102.0/interproscan-data-102.0.tar.gz
tar -pxzf interproscan-data-102.0.tar.gz
```

2. **Run `InterProScan`**

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile <ContainerRuntime,Executor> \
  --input <Fasta> \
  --datadir <Data>
```

For example, to run `InterProScan` locally using Docker:
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile docker,local \
  --input my-fasta.faa \
  --datadir interproscan-data-102.0
```

3. **(Optional) Install licensed software**

By default `Phobius`, `SignalP`, and `TMHMM` member database analyses are deactivated in `InterProScan6` because they contain licensed components. In order to activate these analyses please see the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-tmhmm) documentation.

## Installing from source

1. **Download InterPro data files**

Run the following commands within the `InterProScan` project directory to download and extract all the required data files:

```bash
curl -OJ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/102.0/interproscan-data-102.0.tar.gz
tar -pxzf interproscan-data-102.0.tar.gz
```

2. **Clone the `InterProScan` repository**

```bash
git clone https://github.com/ebi-pf-team/interproscan6.git
```

3. **Pull or build the container images**

Pull down the `InterProScan` image from DockerHub using your container runtime of choice. E.g. using Docker:
```bash
docker pull interpro/interproscan6:latest
```

Or build the Docker image locally:
```bash
docker build -t interproscan6 . <-----
```

Then, optionally, convert the Docker images to alternative container runtimes. For example to convert the Docker image to an Apptainer image:
```bash
docker save interproscan6 > interproscan6.tar
apptainer build interproscan6.sif docker-daemon://interproscan6:latest
```

4. **Run `InterProScan`**

```bash
nextflow run <Path-to-Interproscan-repo>/main.nf \
  -profile <ContainerRuntime,Executor> \
  --input <Fasta> \
  --datadir <Data>
```

For example, to run `InterProScan` locally using Docker:
```bash
nextflow run interproscan6/main.nf \
  -profile docker,local \
  --input my-fasta.faa \
  --datadir interproscan-data-102.0
```

# Using `InterProScan6`

## Quick start

`InterProScan6` is configured via the command-line. The only mandatory arguments are the:

* Runtime profiles (`-profiles`)
* Path to input FASTA file (`--input`)
* Path to the InterPro data (`--datadir`)

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile <container runtime, and executor> \
  --input <path to fasta file> \
  --datadir <path to the data dir>
```

**`-profile` must** be included, and is used to define the executor and the container runtime used.

**`--input` must** be included, and is used to define the path to the input file containing the query sequences to be analysed in FASTA format.

**`--datadir` must** be included, and is used to define the path to the directory containing the downloaded InterPro data.

> [!NOTE]
> The `--datadir` flag is not needed when only running member databases that do not require additional data files.
> This only applies to `mobidblite` and `coils`.

> [!NOTE]  
> Note that `-profile` has a single dash because this is a Nextflow argument, but `--input` and `--datadir` have two dashes because these are `InterProScan6` argument.

## Profiles

A [Nextflow profile](https://nextflow.io/docs/latest/config.html#config-profiles) is a configuration file that defines runtime configuration attributes.

For `InterProScan6` to run, a profile for the container runtime and a profile for an executor must also be defined.

`InterProScan6` includes the following profiles:

* Executor profiles:
  * local
  * slurm
  * lsf
* Container runtime profiles:
  * docker
  * singularity
  * apptainer

For example, to run `InterProScan` on a cluster with the SLURM scheduler and Singularity:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile slurm,singularity \
  --input <path to fasta file> \
  --datadir <data>
```

Nextflow also supports using Charliecloud, Podman, Sarus, and Shifter. However, you will need to build your own [Nextflow profiles](https://nextflow.io/docs/latest/config.html#config-profiles) to support these container runtimes.

## Optional arguments

**`--applications`** - Applications/member databases to run. By default `InterProScan` runs all member databases in the consortium. Use the `--applications` to define a comma separate list of applications names (case insensitive).

**`--disable_precalc`** - `InterProScan6` will check against a remote database of precalculated matches. The downstream analyses will then only be run against query sequences for whom no precalcualted match is available. You can disable this operation and run the analyses against all sequences in the input FASTA file by including the `--disable_precalc` flag in the `InterProScan6` command.

**`--formats`** - List output file formats as a comma separated list. Supported: JSON, TSV, XML. Default: TSV, JSON and XML

**`--goterms`** - Include Gene Ontology (GO) annotations in the final results.

**`--help`** - Display the help message.

**`--max-workers`** - Maximum number of workers available for the `InterProScan` when running locally.
> [!IMPORTANT] 
> - *--max-workers* is only applies when using the `local` profile (i.e. `-profile local`) it does **_not_** apply when running on a cluster.
> - IPS6 will always use a minimum or 2 CPUs, with at least 1 dedicated to the main workflow and 1 to run processes (exception for PRINTS member, which require 2 CPUs to run processes).

**`--nucleic`** - Instead of providing protein sequences, the input FASTA file contains nucleic sequences. See the '[Using DNA sequences](#using-dna-sequences)' section for more information.

**`--outdir`** - Output directory. By default, the output files are written to the current working directory. Use `--outdir` to define the relative or abosolute path to the output directory. `InterProScan` will build all necessary directories.

**`--pathways`** - Include corresponding Pathway annotations in the final results. 

> [!WARNING]
> Nextflow does not tolerate spaces in file paths.

For example, to run `InterProScan6` locally using Docker and only member databases AntiFam and SFLD, without checking for pre-calculated matches in InterPro, and writing the results only in the JSON and XML formats to the directory `results` that include GO terms and pathway annotations (presuming the InterPro data is located in `./data`):

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile docker,local \
  --input files_test/best_to_test.fasta \
  --datadir data \
  --applications signalp,antifam \
  --disable_precalc \
  --formats json,xml \
  --outdir results \
  --goterms \
  --pathways
```

## Applications

Below is a list of the applications (built in and those that require additional installation steps) that are available in `InterProScan6`:

* InterPro Consortium:
  * [Cath-Gene3D]( https://www.cathdb.info/) (use as 'Gene3D' to run Cath-Gene3D in `InterProScan`)
  * [CDD](https://www.ncbi.nlm.nih.gov/cdd)
  * [HAMAP](https://hamap.expasy.org/)
  * [MobiDB-lite](http://old.protein.bio.unipd.it/mobidblite/)
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
  * [Antifam](https://github.com/ebi-pf-team/antifam)
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
> Specifying the eukarya method of `SignalP6` (`SignalP_EUK`) triggers post-processing of the SP predictions by `SignalP6` to prevent spurious results (it will only predict type Sec/SPI).

## Using DNA sequences

`InterProScan6` takes advantage of the Open Reading Frame (ORF) prediction tool `esl-translate` within the [`easel` tool suite](https://github.com/EddyRivasLab/easel) in order to analyse nucleic acid sequences.

The `easel` application itself and all of its dependencies are already integrated in `InterProScan` so no additional set up is required!

**To run anlyses with nucleic acid sequences, run `InterProScan6` with the `--nucleic` flag**

    nextflow run ebi-pf-team/interproscan6 \
        -profile <executor,container runtime> \
        --input <path to fasta file> \
        --data <data-dir> \
        --nucleic

By default `InterProScan6` will assume the input FASTA file contains protein sequences. The `--nucleic` flag instructs `InterProScan6` to retrieve possible ORFs using the `easel` tool suite.

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
> At the moment, only Gene3D and FunFam are able to detect discontinious domains.

### The envelope

The envelope represents the region of a protein sequence where the domain may be located. Often it is wider than what HMMER chooses as a reasonably confident alignment.

### The PIRSF exception

PIRSF uses the envelope as for start and end of the location, and lists the hmmStart and hmmEnd as the aliFrom and aliTo from HMMER3. 

### The Panther exception

The output from HMMER3 against the HMM models of Panther is post-processed to select only the best homologous family. Therefore, there is a maximum of one domain hit for each Panther signature in a protein sequence. Owing to this the E-value and Score are listed under the `signature` key, not the `locations` key.

# Installing licensed applications (`Phobius`, `SignalP`, `TMHMM`)

By default `Phobius`, `SignalP`, and `DeepTMHMM` member database analyses are deactivated in `InterProScan6` because they contain licensed components. In order to activate these analyses please obtain the relevant licenses and files from the provider (ensuring the software version numbers are the same as those supported by your current `InterProScan6` installation).

Files can be placed in any location.

> [!NOTE]  
> The instructions below presume `InterProScan6` is being run with the Docker. If an alternative 
> container runtime is used the methods below will need to be adapated accordingly.
> As above, if using Singularity and Apptainer, the images should be kept in the root of the 
> `InterProScan6` directory. Otherwise, please update the container paths in the respective `utilities/profiles/` config files.

## DeepTMHMM

1. Contact the DeepTMHMM help desk and request a stand-alone licensed copy of DeepTMHMM, and then download the DeepTMHMM `zip` file provided by DeepTMHMM.
2. Unpack the `zip` file
```bash
unzip DeepTMHMM-Academic-License-v1.0.zip -d <TMHMM-DIR>
```
3. Update the TMHMM dir path the `conf/applications.config` file
```groovy
tmhmm {
    name = "tmhmm"
    dir = ""   <--- update the dir path
}
```
4. Run `InterProScan6` with `(Deep)TMHMM`:

With the `dir` field populated in `conf/applications.config`, `TMHMM` will be included as a default application
when the `--application` flag is not used, else include `tmhmm` in the list of applications defined using the `--applications` flag.

For example, when running `InterProScan` locally, using Docker and with the InterPro data dir located at `./data`:
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile local,docker
    --input tests/data/test_prof.fa \
    --datadir data \
    --applications tmhmm
```

## `Phobius`

1. Download Phobius from the [Phobius server](https://software.sbc.su.se/phobius.html)
2. Unpack the `tar` file
```bash
tar -xzf phobius101_linux.tgz -C <PHOBIUS-DIR>
```
3. Build a Docker image using the Dockerfile provided at `utilities/docker/tmhmm/Dockerfile` in the `InterProScan6` repo
```bash
docker image build -t tmhmm utilites/docker/tmhmm
```
4. Update the `Phobius` dir path in the `conf/applications.config` file
```groovy
phobius {
    name = "Phobius"
    invalid_chars = "-*.OXUZJ"
    dir = ""   <--- update the dir path
}
```
5. Run `InterProScan6` with `Phobius`:

With the `dir` field populated in `conf/applications.config`, `Phobius` will be included as a default application 
when the `--application` flag is not used, else include `phobius` in the list of applications defined using the `--applications` flag.

For example, when running `InterProScan` locally, using Docker and with the InterPro data dir located at `./data`:
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile local,docker
    --input tests/data/test_prof.fa \
    --datadir data \
    --applications phobius
```

## `SignalP`

1. Obtain a license and download `SignalP6` (`SignalP` version 6) from the [SignalP6 server](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (under 'Downloads').
    * Either fast or slow models can be implemented
    * To change the implemented mode please see the [Changing mode](#changing-mode) documentation
2. Unpackage the `SignalP6` `tar` file
```bash
tar -xzf signalp-6.0h.fast.tar.gz -C <SIGNALP-DIR>
```
4. Update the dir path in `conf/applications.config`
```groovy
signalp_euk {
    name = "SignalP-Euk"
    organism = "eukarya"
    dir = ""    <-------- update with <SIGNALP-DIR>
    mode = "fast"
}
signalp_prok {
    name = "SignalP-Prok"
    organism = "other"
    dir = ""    <-------- update with <SIGNALP-DIR>
    mode = "fast"
}
```
5. Running `InterProScan6` with `SignalP6`

With the `dir` field populated in `conf/applications.config`, `SignalP` will be included as a default application
when the `--application` flag is not used, else include `signalp_prok` (to set the `SignalP` `--organism` 
argument to `"other"`) and `signalp_euk` (to set the `SignalP` `--organism` argument to `"eukaryote"`) 
in the list of applications defined using the `--applications` flag.

For example, when running `InterProScan` locally, using Docker and with the InterPro data dir located at `./data`:
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile local,docker
    --input tests/data/test_prof.fa \
    --datadir data \
    --applications signalp_euk,signalp_prok
```

> [!NOTE]
> Quoting the SignalP documentation:  
> Specifying the eukarya method of `SignalP6` (`SignalP_EUK`) triggers post-processing of the SP predictions by `SignalP6` to prevent spurious results (it will only predict type Sec/SPI).

### Changing the mode of `Signalp6` in `InterProScan6`

`SignalP6` supports 3 modes: `fast`, `slow` and `slow-sequential`. The mode can be set using the `--signalpMode` flag. The default mode is `fast`.

> [!WARNING]
> The slow mode can take 6x longer to compute. Use when accurate region borders are needed.

You may need to install the other models manually, please see the [SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#installing-additional-modes).

For example, to run `InterProScan` with the input file `best_to_test.fasta`, using SignalP with only eukaryotic models in slow mode, and with retrieving precalculated matches disabled on a local machine using docker, and the InterPro data dir is located at `./data`:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile docker,local \
  --input tests/data/test_prof.fa \
  --datadir data \
  --applications signalp_euk \
  --disable_precalc \
  --signalpMode slow
```

### Run SignalP with GPU acceleration

By default, `SignalP` runs on your CPU. If you have a GPU available, you can convert the `SignalP` models so that 
your installation can use GPU-acceleration. You will need to install `SignalP` in order to convert to GPU models.

1. Convert the ``SignalP`` installation to GPU by following the [SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#converting-to-gpu)
2. Run `InterProScan` with the `--signalpGPU` flag.

For example, to run ``InterProScan`` with only ``SignalP`` enabled, using GPU acceleration on a SLURM cluster with Singularity support:
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile singularity,slurm
  --input <fasta file> \
  --applications signalp \
  --signalpGPU \
```

# Benchmarking and troubleshooting the performance

Nextflow provides some built in options for assessing the operation of `IPS6`, including generating an HTML report. However, these reports are limited to presenting the resource usage from only a single run, and can only be generated if a run is successful. Consequently, we have packaged a simple benchmarking script into IPS6 to enable assessing the task duration and resource usage across multiple runs, and customised grouping of the data. For example, you may wish to clearly see differences in performance with altering the batch size, the number of CPUs or amount of memory allocated. 

You can find the complete details for benchmarking and assessing the performance of `IPS6` in `./benchmarking/README.md`.

In brief:

1. Install all necessary third party packages listed in `benchmarking/requirements.txt`.
2. Build a trace file as part of your InterProScan runs.
3. Update or create the benchmarking JSON config file
4. Run the benchmarking

```bash
# running from the root of the IPS6 project dir
# and using the benchmarking/tracefiles.json file
python3 benchmarking/benchmark_ips6.py benchmarking/tracefiles.json
```

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

# Troubleshooting

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
sudo nextflow run ebi-pf-team/interproscan6 --input <path to fasta file> 
```

Also try providing root privileges to docker within Nextflow, by changing the `runOptions` key in `./utilities/profiles/docker.config` (you will need to be working with a [local installation of `InterProScan`](#installing-from-source) for this fix):

```
 docker {
     enabled = true
     mountFlags = 'Z'
+    runOptions = '--user root'  <<-- ADD
 }
```

## File not found

If you receive a file not found error:

```bash
FileNotFoundError: [Errno 2] No such file or directory
```

This may be due to the docker image not being built correctly. This can happen due to file permissions from preventing docker from reading the `InterProScan6` files.

Try running docker with root privileges:

```bash
sudo docker build -t interproscan6 .
```

Check the docker installation is configured correctly, with all necessary privileges. [StackOverflow](https://stackoverflow.com/questions/48957195/how-to-fix-docker-got-permission-denied-issue)

For example, check root privileges have been provided to the docker socket, although be careful of the security implications of this:

```bash
sudo chmod 666 /var/run/docker.sock
```

## Segmentation fault

If you a receive an error such as the following:

```bash
Command error:
  .command.sh: line 2:     7 Segmentation fault      (core dumped) /opt/hmmer3/bin/hmmsearch --cut_ga --cpu 1 -o 7.0._.antifam._.out AntiFam.hmm mini_test.1.fasta
```

This is generally due to HMMER being unable to find a necessary data file.
Make sure the data directory is correctly structured and populated and `InterProScan` is 
pointed to the correct data directory using the `--data` flag if not using the default data
directory location in the project dir.

## Segmentation fault

If you a receive an error such as the following:
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

A potential fix is to provide root privileges to the docker contains run by Nextflow in `./utilities/profiles/docker.config` (you will need to be working with a [local installation of `InterProScan`](#installing-from-source) for this fix):

```groovy
docker {
    enabled = true
    mountFlags = 'Z'
    runOptions = '--user root'  <---- add this line
}
```
