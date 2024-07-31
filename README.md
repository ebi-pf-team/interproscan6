# InterProScan6

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test) 
![Unit tests](https://github.com/ebi-pf-team/interproscan6/actions/workflows/main.yml/badge.svg)  
[![codecov](https://codecov.io/gh/ebi-pf-team/interproscan6/graph/badge.svg?token=7MP9WCJHAQ)](https://codecov.io/gh/ebi-pf-team/interproscan6)
[![Citation](https://github.com/ebi-pf-team/interproscan6/actions/workflows/citation.yml/badge.svg)](#citation)

**!! UNDER DEVELOPMENT !!**

[InterPro](http://www.ebi.ac.uk/interpro/) is a database which integrates together predictive information about proteins’ function from a number of partner resources, giving an overview of the families that a protein belongs to and the domains and sites it contains.

Users who have novel nucleotide or protein sequences that they wish to functionally characterise can use the software package `InterProScan` to run the scanning algorithms from the InterPro database in an integrated way. Sequences are submitted in FASTA format. Matches are then calculated against all of the required member database’s signatures and the results are then output in a variety of formats.

## Table of Contents
<!-- TOC -->
- [`InterProScan`](#interproscan6)
- [Setup and Requirements](#set-up-and-requirements)
  - [Requirements](#requirements)
  - [Set up](#set-up)
- [Using `InterProScan6`](#using-interproscan6)
  - [Quick start](#quick-start)
  - [Using DNA sequences](#using-dna-sequences)
  - [Inputs and configuration](#inputs-parameters-and-configuration)
  - [Outputs and results](#outputs)
- [Installing licensed applications (`Phobius`, `SignalP`, `TMHMM`)](#installing-licensed-applications-phobius-signalp-tmhmm)
- [Citing `InterProScan`](#citation)
- [Trouble shooting](#trouble-shooting)
<!-- /TOC -->

# Set up and requirements

## Requirements

* `Java` (version >= 11)
* `Nextflow` (version >=23.04.02)
* `Docker` (version >= 24.0.5)

## General set up

1. Download member data files:

    curl ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.67-99.0/alt/interproscan-data-5.67-99.0.tar.gz \
        --output interproscan-data-5.67-99.0.tar.gz
    tar -pxzf interproscan-data-5.67-99.0.tar.gz
    mv interproscan-5.67-99.0/data .
    rm interproscan-5.67-99.0 -rf
    rm interproscan-data-5.67-99.0.tar.gz

These commands download and store all member database data in the `data/` directory.

If you do not store these data in another directory you will need to update the members config file: 
`subworkflows/sequence_analysis/members.config`. Specifically, you will need to 
update the paths for:
* `hmm`
* Under `postprocess`: [_where applicable_]
  * `data`
  * `models`
  * `modesl_dir`
  * `rules`

Please provide absolute paths. You can use the `$projectDir` short cut to represent the path to the root directory of `InterProScan6`.

2. Download InterPro data files (xref (entries, goterms and pathways)):

    mkdir i6data
    python interproscan6/data_xrefs/get_data_to_i6.py  # for devs

This stores the InterPro data in the `data_xrefs/` directory.

If you store these data in an alternative directory, please update the paths under the `xrefs` 
heading in the corresponding configuration file (`./interproscan/subworkflows/xrefs/xrefs.config`).

Please provide absolute paths. You can use the `$projectDir` short cut to represent the path to 
the root directory of `InterProScan6`.

3. Build a docker `InterProScan6` base image (this includes all non-licensed dependencies including `HMMER`, `BLAST`, `BioPython`, `easel`, etc.)

    docker build -t interproscan6 .

4. [Optional] install licensed software

By default `MobiDB`, `Phobius`, `SignalP`, and `TMHMM` member database analyses are deactivated in `InterProScan6` 
because they contain licensed components. In order to activate these analyses please see the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-tmhmm) documentation.

## Singularity set up

Not all systems support using Docker, therefore, the `interproscan6` Docker image will need to be 
converted to another virtualization system. Singularity is an alternative container runtime to
Docker that does not require root privileges or a separate daemon process, unlike Docker.

You will need Singularity installed on the system you will use to run the pipeline.

1. Follow the general set up laid out above on the system you are going to run the pipeline. Build
the `interproscan6` Docker image (in step 3 of the general set up) on a system with Docker enabled -
you **do not** need have downloaded the InterPro release data in order to build the Docker image.

2. Save the `interproscan6` Docker image to a `tar` archive.

```bash
docker save interproscan6 > interproscan6.tar
```

3. Build a Singularity image from the archived Docker image:

```bash
singularity build interproscan6.img docker-archive://interproscan6.tar
```

4. Test the image on the system you will use to run the pipeline.

```bash
singularity shell interproscan6.img
  $ ls
  $ exit
```

Keep the Singularity image (`interproscan6.img`) in the root of the `InterProScan6` repository 
on the sytem you will use to run the pipeline. For example, you can create an `interproscan6.img` on your 
local laptop, and upload it to the HPC you will use to run `InterProScan6`.

When running `InterProScan6` with Singularity include `singularity` in the `-profiles` option. E.g.:
```bash
nextflow run interproscan.nf --input my_seqs.fasta -profiles singularity
```
The order the profiles are listed after `-profiles` does **not** matter.

# Using `InterProScan6`

## Quick start

How to run:

    nextflow run interproscan.nf --input <path to fasta file> -profile <container runtime, and executor>

If running `InterProScan6` locally, you need only provide the corresponding container runtime 
profile. Currently, Docker (profile: `docker`), Singularity (profile: `singularity`), and Apptainer 
(profile: `apptainer`) are supported.

If running `InterProScan6` on a cluster you will need to *additionally* supply the executor profile,
for example, for running on a cluster using the SLURM scheduler, add the SLURM profile, for example:

```bash
  nextflow run interproscan.nf 
    --input <path to fasta file> 
    -profile slurm,singularity
```

The results will apear in `result/` folder.  
For debugging, you an find all working files generated by the pipeline in the `work/` dir.

An example command to run `InterProScan6`, only using the `AntiFam` member database and 
`SignalP`, without checking for pre-calculated matches in InterPro (using an example input file), 
using a local system and Docker:

```bash
  nextflow run interproscan.nf \
    --input files_test/best_to_test.fasta \
    --applications signalp,antifam \
    --disable_precalc \
    -profile docker
```

## Using DNA sequences

`InterProScan6` takes advantage of the Open Reading Frame (ORF) prediction tool `esl-translate` within the [`easel` tool suite](https://github.com/EddyRivasLab/easel).

The `easel` application itself and all of its dependencies are integrated in InterProScan via a docker image (see [set-up](#set-up)).

### Quick start

Run `InterProScan6` with the `--nucleic` flag

    nextflow run interproscan.nf \
        --input <path to fasta file> \
        --nucleic

By default `InterProScan6` will assume the input FASTA file contains protein sequences. The `--nucleic` flag instructs `InterProScan6` to retrieve all possible ORFs using the `easel` tool suite.

**Note:** The input FASTA file must contain sequences of the same type, i.e. _all_ protein sequences or _all_ nucleic acid sequences.

### Configure

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

## Inputs, parameters and configuration

`InterProScan6` is configured via the command-line. The only mandatory parameter is `--input`.

**IMPORTANT:** For parameters that have more than one value, separate values using a comma (not spaces) (e.g. `--applications antifam,ncbifam,pfam)`)

**Configuration parameters:**

* `--applications` - a list of member databases/applications to be employed in the analysis. By default, `InterProScan6` employs all member databases and applications, use this flag to run only a subset of applications. For example:

    nextflow run interproscan.nf --input <fasta file> --applications NCBIfam,Panther,Pfam

```yaml
applications: AntiFam,CDD,Coils,FunFam,Gene3d,HAMAP,MobiDBLite,NCBIfam,Panther,Pfam,Phobius,PIRSF,PIRSR,PRINTS,PrositePatterns,PrositeProfiles,SFLD,SignalP_EUK,SignalP_GRAM_NEGATIVE,SignalP_GRAM_POSITIVE,SMART,SuperFamily,TMHMM
```

* `--disable_precalc` - Do not run comparison against an InterPro release to retrive precalculated matches, instead run `interproscan` for all input sequences. [Boolean]
* `--formats` - List output file formats. Supported: JSON,TSV,GFF,XML
* `--goterms` - Whether to retrieve and include Gene Ontology terms from InterPro in the output files. [Boolean]
* `--help` - Whether to disble the help message - `InterProScan6` will not run any analysis when `help` is set to true. [Boolean]
* `--input` - Path to input FASTA file
* `--outfile` - Path and prefix for output files
* `--pathways` - Optional, switch on lookup of corresponding Pathway annotation (IMPLIES - `lookup_file` is defined) [Boolean]
* `--lookup_file` - Lookup of corresponding InterPro annotation in the TSV and GFF3 output formats.

### Nextflow configuration

Configure the `InterProScan6` utility operations by updating `./nextflow.config`:

* Define the batch size (number of sequences included in each batch job)
* Change the accepted member databases (important if including `SignalP`, `TMHMM` and `Phobius`)
* Update the paths to the InterPro release files

### Sequences

The input sequences to be analysed by `InterProScan6`must be provided in FASTA format. All input sequences must be provided in a multi-sequence FASTA file.

At the moment only protein (amino acid) sequences are supported.

## Outputs

:TODO:

# Installing licensed applications (`Phobius`, `SignalP`, `TMHMM`)

By default `Phobius`, `SignalP`, and `TMHMM` member database analyses are deactivated in `InterProScan6` 
because they contain licensed components. In order to activate these analyses please 
obtain the relevant license and files from the provider (ensuring the software version 
numbers are the same as those supported by your current `InterProScan6` installation).

Files can be placed in any location.

## `MobiDB`

Some of the compoments within `MobiDBLite` are GPL-licensed, meaning all software and data, and thus 
work that uses this software, also needs to be GPL-licensed. This may not be ideal or suitable
for all users. Therefore, we provide a version of the `MobiDBLite` analytical software that 
is not GPL-licensed, called `idpred`.

To setup `MobiDB`/`idpred` for InterProScan6 you only need to create a Docker image using the 
provided Dockerfile: From the root of this repository:

```bash
cd docker_files/mobidb
docker build -t idrpred .
```

## `SignalP`

### Adding `SignalP` (version 6) to `InterProScan6`

1. Obtain a license and download `SignalP6` (`SignalP` version 6) from the `SignalP6` [server](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (under 'Downloads').
    * Either fast or slow models can be implemented
    * To change the implemented mode please see the [Changing mode](#changing-mode) documentation

2. Unpackage the `SignalP6` `tar` file

    tar -xzf signalp-6.0h.fast.tar.gz -C <SIGNALP-DIR>

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

3. Update the `InterProScan6` configuration: specifically, update `subworkflows/sequence_analysis/members.config`:
```
signalp {
    release = "6.0h"  <--- make sure the release is correct
    runner = "signalp"
    data {
        mode = "fast"    <--- UPDATE MODE: fast, slow, or slow-sequential
        model_dir = "$projectDir/bin/signalp/models"  <--- UPDATE PATH TO models DIR
        organism = "other"
    }
}
```

**Note:** _Set `organism` in `nextflow.config` to `"eukaryote"` or `"euk"` to limit the predictions to Sec/SPI, or leave as `"other"` to apply all models in `SignalP6`, as per the `SignalP6` documentation:_

> Specifying the eukarya method of `SignalP6` (`SignalP_EUK`) triggers post-processing of the SP predictions by `SignalP6` to prevent spurious results (only predicts type Sec/SPI).

4. Add `SignalP` to the application list in `nextflow.config`:

params {
    batchsize = 100
    help = false
    applications = 'AntiFam,CDD,Coils,FunFam,Gene3d,HAMAP,MobiDBLite,NCBIfam,Panther,Pfam,PIRSF,PIRSR,PRINTS,PrositePatterns,PrositeProfiles,SFLD,SMART,SuperFamily,SignalP' <--- ADD NEW APPLICATION
    disable_precalc = false
}

### Running `InterProScan6` with `SignalP6` enabled

Include `signalp` in the list of applications defined using `--applications` flag.

For example:

    nextflow run interproscan.nf --input files_test/best_to_test.fasta --applications signalp --disable_precalc

### Changing mode

`SignalP6` supports 3 modes: `fast`, `slow` and `slow-sequential`. To change the mode of `SignalP6`:

1. Incorporate the new mode into your `SignalP6` installtion as per the `SignalP6` [documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#installing-additional-modes).

2. Update the `member.config` configuration (`subworkflows/sequence_analysis/members.config`)
```
signalp {
    runner = "signalp"
    data {
        mode = "fast"    <--- UPDATE MODE: fast, slow, or slow-sequential
        model_dir = "$projectDir/bin/signalp/models"
        organism = "other"
    }
}
```

**Note:** _`InterProScan6` only supports the implementation of one `SignalP` mode at a time. A separate `InterProScan6` but be completed for each mode of interest, in order ro apply multiple modes to the same dataset_.

### Converting from CPU to GPU, and back again

:TODO: -- add support for GPU run

The model weights that come with the `SignalP` installation by default run on your CPU.
If you have a GPU available, you can convert your installation to use the GPU instead. 

1. Convert the `SignalP` installation to GPU by following the `SignalP` [documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#converting-to-gpu)
2. ....???....

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
```
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

    sudo docker build -t interproscan6 .


Check the docker installtion is configured correctly, with all necessary privileges. [StackOverflow](https://stackoverflow.com/questions/48957195/how-to-fix-docker-got-permission-denied-issue)

For example, check root privileges have been provided to the docker socket

    sudo chmod 666 /var/run/docker.sock


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
