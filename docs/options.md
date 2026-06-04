<style>
    .md-typeset__table {
        min-width: 100%;
    }

    .md-typeset table:not([class]) {
        display: table;
    }

    table td:first-child,
    table th:first-child {
        white-space: nowrap;
    }
</style>

# Options / Parameters

InterProScan 6 is distributed as a Nextflow workflow, so a typical command combines:

- Nextflow options with a single dash, such as `-profile`, `-r`, and `-w`
- InterProScan workflow parameters with a double dash, such as `--input` and `--formats`

The tables below describe the InterProScan 6 parameters, plus the main Nextflow options that are part of almost every invocation.

## Core invocation options

| Option | Description |
|--------|-------------|
| `-r`, `-revision <VERSION>` | Tells Nextflow which InterProScan version to run. Recommended for reproducibility. For example, use `-r 6.0.1` to pin a specific InterProScan release. |
| `-profile <PROFILE>` | Selects the runtime environment, such as `docker`, `slurm`, etc. Common combinations include `docker,test` and `singularity,slurm`. |
| `--input <FASTA>` | Path to the input FASTA file. Required unless a profile such as `test` supplies it automatically. |
| `--datadir <DATADIR>` | Directory used for InterPro and member-database data files. |
| `--interpro <VERSION>` | InterPro release to use. The default is `latest`. A pinned version such as `108.0` is better for reproducibility. |

!!! tip

    `--datadir` is only required when at least one selected analysis needs external data files. It is not required if you run only self-contained analyses such as COILS, MobiDB-lite, or TMbed.

## Input interpretation

| Option | Description |
|--------|-------------|
| `--nucleic` | Treats the input as nucleotide FASTA. The workflow translates all six reading frames and analyses the resulting ORFs. Disabled by default; protein FASTA is assumed otherwise. |

## Analysis selection

| Option | Description |
|--------|-------------|
| `--applications <LIST>` | Comma-separated list of analyses to run. By default, InterProScan runs all analyses except the deep-learning-based analyses. |
| `--skip-applications <LIST>` | Comma-separated list of analyses to exclude. The default is to skip none. |
| `--run-ml` | Enables deep-learning-based analyses such as DeepTMHMM, SignalP, and TMbed. Disabled by default because these analyses are more resource-intensive. |

!!! info

    Some analyses rely on separately licensed third-party applications. This affects the external software component, not your ability to run InterProScan on your own sequences. See [Licensed applications](licensed_applications.md).

!!! tip

    Analysis names passed to `--applications` and `--skip-applications` are case-insensitive, and hyphens, underscores, and spaces are ignored.

!!! warning

    - `--applications` cannot be combined with `--skip-applications`.
    - `--applications` cannot be combined with `--run-ml`.

## Output control

| Option | Description |
|--------|-------------|
| `--outdir <DIR>` | Directory where output files are written. The default is the current working directory (`.`). |
| `--outprefix <PREFIX>` | Base name for output files, without directory or extension. The default is the input filename. It must not contain slashes. |
| `--formats <LIST>` | Comma-separated output formats. Supported values are `gff3`, `json`, `jsonl`, `tsv`, and `xml`. The default is all formats. |
| `--goterms` | Includes Gene Ontology annotations in supported output files. Disabled by default. |
| `--pathways` | Includes pathway annotations in supported output files. Disabled by default. |
| `--skip-repr-locations` | Skips the representative-location selection step. Disabled by default. |

!!! info

    `--outprefix` controls only the filename stem. Use `--outdir` to choose the output directory.

## Execution and resource control

| Option | Description |
|--------|-------------|
| `--cpus <N>` | Number of CPU threads allocated per task for applications that support multithreading. The default is `1`. |
| `--use-gpu` | Enables GPU acceleration for compatible deep-learning-based analyses. This only has an effect when GPU-capable analyses are selected and configured. |
| `--max-workers <N>` | Maximum number of parallel tasks when running locally. If unset, Nextflow uses its normal scheduler behaviour. |

## External services and data

| Option | Description |
|--------|-------------|
| `--matches-api-url <URL>` | Uses an alternative InterPro Matches API endpoint for precomputed results. Default: [https://www.ebi.ac.uk/interpro/matches/api](https://www.ebi.ac.uk/interpro/matches/api). |
| `--no-matches-api` | Disables the Matches API and forces all analyses to run locally. Disabled by default. |
| `--globus` | Downloads InterProScan data from the Globus mirror instead of the EMBL-EBI FTP site. Disabled by default. |
| `--skip-interpro-version-check` | Skips the compatibility check between the workflow version and the selected InterPro data release. Disabled by default. |

## Other options

| Option | Description |
|--------|-------------|
| `--help` | Prints the usage message and exits |
