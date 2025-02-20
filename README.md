# InterProScan6

[![nf-test](https://img.shields.io/badge/tested_with-nf--test-337ab7.svg)](https://github.com/askimed/nf-test)
![Unit tests](https://github.com/ebi-pf-team/interproscan6/actions/workflows/unit-tests.yml/badge.svg)
[![Check Docker Files](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml/badge.svg)](https://github.com/ebi-pf-team/interproscan6/actions/workflows/docker-check.yml)
[![Citation](https://github.com/ebi-pf-team/interproscan6/actions/workflows/citation.yml/badge.svg)](#citation)

> [!CAUTION]
> InterProScan6 is currently under active development and is not yet stable enough for a full release.

[InterPro](http://www.ebi.ac.uk/interpro/) is a database which integrates together predictive information about proteins’ functions from a number 
of partner resources, giving an overview of the families that a protein belongs to as well as the domains and sites it contains.

Users who have novel nucleotide or protein sequences that they wish to functionally characterise can use the software 
package `InterProScan` to run the scanning algorithms from the InterPro database in an integrated way. Sequences are 
submitted in FASTA format. Matches are then calculated against all the required member databases signatures and the 
results are then output in a variety of formats.

# Documentation

Our full documentation is available at [ReadTheDocs](https://interproscan-docs.readthedocs.io/en/latest/).

# Requirements

* `Nextflow` (version >=23.04.01)
* A container run time. `InterProScan` includes built in support for:
    * `Docker` (version >= 24.0.5)
    * `SingularityCE` (version >= 4.2.0)
    * `Apptainer` (version >= 1.3.4)
* Licenses and additional data from their respective authors are required to run `Phobius`, `SignalP` and `DeepTMHMM`

# Set up

## Quick

Nextflow pulls the latest `InterProScan` image from DockerHub if one is not already available on the host system,
otherwise Nextflow will use the available image.

1. **Download InterPro data files**
```bash
curl -OJ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0/104.0/interproscan-data-104.0.tar.gz
tar -pxzf interproscan-data-104.0.tar.gz
```

2. (Optional) Install licensed software (`Phobius`, `SignalP`, `DeepTMHMM`) - See the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-deeptmhmm) documentation.

3. **Run `InterProScan`**
```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile <ContainerRuntime (docker,singularity,apptainer), Executor (local,slurm,lsf)> \
  --input <Fasta> \
  --datadir <Data>
```

## Installing from source

1. **Download InterPro data files**
```bash
curl -OJ https://ftp.ebi.ac.uk/pub/software/unix/iprscan/6/6.0/104.0/interproscan-data-104.0.tar.gz
tar -pxzf interproscan-data-104.0.tar.gz
```

2. **Clone the `InterProScan` repository**
```bash
git clone https://github.com/ebi-pf-team/interproscan6.git
```

3. **Pull or build the container images**

Pull from DockerHub using a container runtime. E.g. using Docker:
```bash
docker pull interpro/interproscan6:latest
```
Or build the Docker image locally (then optionally convert to an alternative container runtime):
```bash
docker build -t interproscan6 utilities/docker/interproscan
```

4. (Optional) Install licensed software (`Phobius`, `SignalP`, `DeepTMHMM`) - See the ['Installing licensed applications'](#installing-licensed-applications-phobius-signalp-deeptmhmm) documentation.

5. **Run `InterProScan`**
```bash
nextflow run <Path to the InterProScan main.nf file> \
  -profile <ContainerRuntime (docker,singularity,apptainer), Executor (local,slurm,lsf)> \
  --input <Fasta> \
  --datadir <Data>
```

# How to run `InterProScan6`

`InterProScan6` is configured via the command-line. The only mandatory arguments are the:

* Runtime profiles (`-profiles`) (an executor (local, slurm or lsf) and a container (docker, singularity or apptainer))
* Path to input FASTA file (`--input`)
* Path to the InterPro data (`--datadir`)

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile <docker,singularity,apptainer...local,slurm,lsf> \
  --input <path to fasta file> \
  --datadir <path to the data dir>
```

> [!WARNING]
> Nextflow does not tolerate spaces in file paths.

> [!NOTE]
> The `--datadir` flag is not needed when only running member databases that do not require additional data files.
> This only applies to `mobidblite` and `coils` (which do not require additional datafiles) and the licensed software
> (`SignalP`, `Phobius`, and `DeepTMHMM`).

**Optional arguments:**

* `--applications` - Comma separated list of analyses to run. By default `InterProScan` runs all members in `conf/applications.conf` with a populated `dir` field.
* `--offline` - Do not retrieve pre-calculated matches from the InterPro Matches API.
* `--formats` - Comma separated list of output files to write (TSV, JSON, XML). Default: TSV, JSON and XML
* `--goterms` - Include Gene Ontology (GO) annotations in the final results.
* `--help` - Display the help message.
* `--max-workers` - Maximum number of workers available for the `InterProScan` when running locally.
* `--nucleic` - The input FASTA file contains nucleic sequences. See the '[Using DNA sequences](#using-dna-sequences)' section for more information.
* `--outdir` - Output directory. By default, the output files are written to the current working directory. `InterProScan` will build all necessary directories.
* `--pathways` - Include corresponding Pathway annotations in the final results.

> [!IMPORTANT]
> *--max-workers* is only applies when using the `local` profile (i.e. `-profile local`) it does **_not_** apply when running on a cluster.
> IPS6 will always use a minimum or 2 CPUs, with at least 1 dedicated to the main workflow and 1 to run 
> processes (exception for PRINTS member, which require 2 CPUs to run processes).

For example, to run `InterProScan6` locally using Docker and only member databases AntiFam and SFLD, without checking 
for pre-calculated matches in InterPro, and writing the results only in the JSON and XML formats to the directory 
`results` that include GO terms and pathway annotations (presuming the InterPro data is located in `./data`):

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile docker,local \
  --input tests/data/test_prot.fa \
  --datadir data \
  --applications signalp,antifam \
  --offline \
  --formats json,xml \
  --outdir results \
  --goterms \
  --pathways
```

## Using DNA sequences

To run anlyses with nucleic acid sequences, run `InterProScan6` with the `--nucleic` flag:

```bash
nextflow run ebi-pf-team/interproscan6 \
    -profile <executor,containerRuntime> \
    --input <path to fasta file> \
    --data <data-dir> \
    --nucleic
```

# Installing licensed applications (`Phobius`, `SignalP`, `DeepTMHMM`)

By default `Phobius`, `SignalP`, and `DeepTMHMM` analyses are deactivated in 
`InterProScan6` because they contain licensed components.

For these software, with their `dir` field populated in `conf/applications.config` they will be included as a 
default application (i.e. when `--application` is not used).

1. Obtain a license and download the analysis software:
   * `DeepTMHMM`: Contact the DeepTMHMM help desk and request a stand-alone licensed copy of DeepTMHMM.
   * `Phobius`: From the [Phobius server](https://software.sbc.su.se/phobius.html
   * `SignalP`: From the [SignalP6 server](https://services.healthtech.dtu.dk/services/SignalP-6.0/) (under 'Downloads')
2. Unpack the `ZIP` or `TAR` file
3. Update the respective members `dir` path in `conf/applications.conf`
```groovy
deeptmhmm {
    name = "DeepTMHMM"
    dir = ""      <--- update the dir path
}
phobius {
  name = "Phobius"
  invalid_chars = "-*.OXUZJ"
  dir = ""        <--- update the dir path
}
signalp_euk {
  name = "SignalP-Euk"
  organism = "eukarya"
  dir = ""        <--- update the dir path
  mode = "fast"
}
signalp_prok {
  name = "SignalP-Prok"
  organism = "other"
  dir = ""        <--- update the dir path
  mode = "fast"
}
```

Quoting the SignalP documentation: "Specifying the eukarya method of `SignalP6` 
(`SignalP_EUK` in `InterProScan6`) triggers post-processing of the SP predictions by 
`SignalP6` to prevent spurious results (it will only predict type Sec/SPI)."

## Changing the mode of `Signalp6` in `InterProScan6`

`SignalP6` supports 3 modes: `fast`, `slow` and `slow-sequential`. The mode can be set using the `--signalpMode` flag. The default mode is `fast`.
You may need to install the other models manually, please see the 
[SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#installing-additional-modes).

```bash
nextflow run ebi-pf-team/interproscan6 \
  -profile <docker,singularity,apptainer...local,slurm,lsf> \
  --input <fasta-path> \
  --datadir <data-dir-path> \
  --signalpMode slow
```

> [!WARNING]
> The slow mode can take 6x longer to compute. Use when accurate region borders are needed.

## Run SignalP with GPU acceleration

1. Convert the ``SignalP`` installation to GPU by following the [SignalP documentation](https://github.com/fteufel/signalp-6.0/blob/main/installation_instructions.md#converting-to-gpu)
2. Run `InterProScan` with the `--signalpGPU` flag.

# Citation

If you use `InterProScan` or `InterPro` in your work, please cite the following publications:

**`InterProScan`:**

> Jones P, Binns D, Chang HY, Fraser M, Li W, McAnulla C, McWilliam H, Maslen J, Mitchell A, Nuka G, Pesseat S, Quinn AF, Sangrador-Vegas A, Scheremetjew M, Yong SY, Lopez R, Hunter S. InterProScan 5: genome-scale protein function classification. Bioinformatics. 2014 May 1;30(9):1236-40. doi: 10.1093/bioinformatics/btu031. Epub 2014 Jan 21. PMID: 24451626; PMCID: PMC3998142.

**InterPro:**

> Paysan-Lafosse T, Blum M, Chuguransky S, Grego T, Pinto BL, Salazar GA, Bileschi ML, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunić I, Marchler-Bauer A, Mi H, Natale DA, Orengo CA, Pandurangan AP, Rivoire C, Sigrist CJA, Sillitoe I, Thanki N, Thomas PD, Tosatto SCE, Wu CH, Bateman A. InterPro in 2022. Nucleic Acids Res. 2023 Jan 6;51(D1):D418-D427. doi: 10.1093/nar/gkac993. PMID: 36350672; PMCID: PMC9825450.

You can find more information about the InterPro consortium (including the member databases) on the [InterPro website](https://www.ebi.ac.uk/interpro/about/consortium/).

# Troubleshooting

For further assistance with installing and using InterProScan, please reach out to us through our 
[help desk](http://www.ebi.ac.uk/support/interproscan) or create an [issue](https://github.com/ebi-pf-team/interproscan6/issues) on our GitHub repository.
