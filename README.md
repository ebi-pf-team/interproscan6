# InterProScan 6

[InterPro](http://www.ebi.ac.uk/interpro/) is a database which integrates together predictive information about proteins’ functions from a number of partner resources, giving an overview of the families that a protein belongs to as well as the domains and sites it contains.

**InterProScan** is the software package that allows sequences to be scanned against InterPro's member database signatures. Researchs who have novel nucleotide or protein sequences that they wish to functionally characterise can use InterProScan to run the scanning algorithms against the InterPro database in an integrated way.

> [!CAUTION]
> InterProScan6 is under active development and is not guaranteed to be stable.

## Installation

To run InterProScan, you need:

* [Nextflow](https://www.nextflow.io/) 24.10.04 or later
* A container runtime. InterProScran currently supports:
    * [Docker](https://www.docker.com/)
    * [SingularityCE](https://sylabs.io/singularity/)
    * [Apptainer](https://apptainer.org/)
* Licenses and additional data from their respective authors are required to run Phobius, SignalP and DeepTMHMM (see [licensed analyses](#licensed-analyses))

## Usage

To test InterProScan, run the following command:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile test,docker \
  --datadir data \
  --interpro latest \
  --download
```

Explanation of parameters:

* `-profile test,docker`: use the following built-in profiles:
  * `test`: use an included FASTA file of protein sequences
  * `docker`: execute tasks in Docker containers
* `--datadir data`: use `data` as the directory to store files required to scan sequences, such as HMM databases. The directory is created if it doesn't exists.
* `--interpro latest`: use the latest version of InterPro and its member databases
* `--download`: download InterPro metadata and member database files, if they are not found in the `data` directory.

InterProScan will create three files in your current working directory:

* `test.faa.json`: sequence annotations in the JSON format
* `test.faa.tsv`: sequence annotations in the TSV format
* `test.faa.xml`: sequence annotations in the XML format

The JSON and XML output files are the richer, while the TSV provides less information.

### Using your own input file

To annotate your own sequences, pass `--input /path/to/your/sequences.fasta` when running InterProScan. Remember to run InterProScan *without* the `test` profile.

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile docker \
  --datadir data \
  --interpro latest \
  --input /path/to/sequences.fasta
```

Got nucleic sequences instead of protein ones? Pass the `--nucleic` parameter.

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile docker \
  --datadir data \
  --interpro latest \
  --input /path/to/sequences.fna \
  --nucleic  
```

### Running specific analyses

One may run only a specific subset of analyses available in InterProScan, e.g. only Pfam for protein family classification and MobiDB-lite for intrinsically disordered regions predicton.

To do so, pass the `--applications` parameter, followed by the comma-separated list of analyses to execute:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -profile test,docker \
  --datadir data \
  --interpro latest \
  --applications Pfam,MobiDB-lite
```

Please note that analyses are case-insensitive and that dashes (`-`) are removed, thereby `MobiDB-lite`  and `mobidblite` are both valid.

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

> [!IMPORTANT]  
> `DeepTMHMM`, `Phobius`, `SignalP-Euk`, and `SignalP-Prok` contain licensed components that prevent them to be integrated by default in InterProScan. Please read [how to execute licensed analyses](#licensed-analyses).

> [!TIP]
> When executing exclusively analyses that do not require InterPro data, the `--datadir`, `--interpro`, and `--download` aren't used. These analyses are `COILS`, `DeepTMHMM`,`MobiDB-lite`, `Phobius`, `SignalP-Euk`, and `SignalP-Prok`.

### Licensed analyses

DeepTMHMM, Phobius, and SignalP include licensed components that cannot be included with InterProScan, therefore these analyses are deactivated by default.

To enable and execute these analyses, you first need to obtain their license, and data.

#### Obtaining licensed components

For each analyses, you need to request a license, then download and extract an archive that contains the data and machine learning model(s) required to run the analysis.

##### DeepTMHMM 1.0

Request a standalone copy of DeepTMHMM 1.0 by sending an email to <licensing@biolib.com>, then extract it when you receive it:

```sh
unzip -q DeepTMHMM-v1.0.zip
```

and make not of the full path to the package:

```sh
echo "${PWD}/DeepTMHMM
```

##### Phobius 1.01

Download a copy of Phobius 1.01 [from Erik Sonnhammer's website](https://software.sbc.su.se/phobius.html), then extract it:

```sh
tar -zxf phobius101_linux.tgz
```

and make not of the full path to the package:

```sh
echo "${PWD}/phobius"
```

##### SignalP 6.0

SignalP 6.0 supports two modes: a slow one that uses the full model, and a fast one that uses a distilled version of the full mode. InterProScan suppports both, but only one at a time. The fast mode is recommended for most users.

You need a license for each of these mode:

* download the distilled model: [fast](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=fast)
* download the full model: [slow_sequential](https://services.healthtech.dtu.dk/cgi-bin/sw_request?software=signalp&version=6.0&packageversion=6.0h&platform=slow_sequential)

Extract the archive:

```sh
tar -zxf signalp-6.0h.fast.tar.gz
```

and make not of the full path to the package:

```sh
echo "${PWD}/signal6p_fast"
```

#### Executing licensed analyses

Once you have downloaded and extracted or or multiple licensed analyses, create a config file as follows:

```groovy
# licensed.conf
params {
    appsConfig {
        deeptmhmm {
            dir = "/path/to/DeepTMHMM"
        }
        phobius {
            dir = "/path/to/phobius"
        }
        signalp_euk {
            dir = "/path/to/signal6p_fast"
        }
        signalp_prok {
            dir = "/path/to/signal6p_fast"
        }
    }
}
```

Then pass this config file with `-c licensed.conf` when running InterProScan:

```sh
nextflow run ebi-pf-team/interproscan6 \
  -c licensed.conf \
  -profile docker
  --input /path/to/sequences.fasta \
  --applications deeptmhmm,phobius,signalp6_euk,signalp6_prok
```

> [!NOTE]  
> SignalP 6.0 has an option to post-process signal peptide predictions to prevent spurious results in eukaryotic proteins. Running InterProScan with `signalp6_euk` *and* `signalp6_prok` will execute SignalP twice, once with the post-processing (better for eukaryotic proteins) and once without (better for prokaryotic proteins). You may want to run only one.

## Documentation

Our full documentation is available at [ReadTheDocs](https://interproscan-docs.readthedocs.io/en/v6/).

# Citation

If you use InterPro in your work, please cite the following publication:

> Blum M, Andreeva A, Florentino LC, Chuguransky SR, Grego T, Hobbs E, Pinto BL, Orr A, Paysan-Lafosse T, Ponamareva I, Salazar GA, Bordin N, Bork P, Bridge A, Colwell L, Gough J, Haft DH, Letunic I, Llinares-López F, Marchler-Bauer A, Meng-Papaxanthos L, Mi H, Natale DA, Orengo CA, Pandurangan AP, Piovesan D, Rivoire C, Sigrist CJA, Thanki N, Thibaud-Nissen F, Thomas PD, Tosatto SCE, Wu CH, Bateman A. **InterPro: the protein sequence classification resource in 2025**. *Nucleic Acids Res*. 2025 Jan;53(D1):D444-D456. [doi: 10.1093/nar/gkae1082](https://doi.org/10.1093/nar/gkae1082).

# Support

For further assistance, please [create an issue](https://github.com/ebi-pf-team/interproscan6/issues) or [contact us](http://www.ebi.ac.uk/support/interproscan).
