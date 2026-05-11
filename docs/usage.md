# Running InterProScan

## Quickstart

The typical command for running InterProScan is as follows:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \ # (1)!
  -profile docker \ # (2)!
  --input /path/to/sequences.faa \ # (3)!
  --datadir data \ # (4)!
  --interpro latest # (5)!
```

1. Run InterProScan `v6.0.0`
2. Executes tasks in Docker containers
3. Path to your input FASTA file of protein sequences
4. Directory where to download and save the reference files
5. Use the latest version of InterPro

!!! info

    Nextflow options use a single dash (e.g. `-profile`, `-r`), while InterProScan parameters use double dashes (e.g. `--input`, `--datadir`).

!!! tip

    `--interpro <VERSION>` is optional and defaults to `latest`, but we strongly recommend to pin a specific version (e.g. `--interpro 108.0`) for reproducibility.

## Selecting specific analyses

By default, InterProScan runs all analyses except the deep-learning-based analyses.

You can control the analyses to be executed with the following parameters:

- `--applications <LIST>`: run only selected analyses
- `--skip-applications <LIST>`: exclude analyses
- `--run-ml`: enable deep-learning-based analyses

`<LIST>` is a comma-separated list of analysis names, e.g. `AntiFam,CATH-Gene3D,Pfam`.

See [Analyses](analyses.md) for the full list of supported analyses.

!!! warning

    `--applications` is mutually exclusive with `--skip-applications` and `--run-ml`.


```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa \
  --applications Pfam,MobiDB-lite # (1)!
```

1. Only Pfam and MobiDB-lite will be executed

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa \
  --skip-applications CDD,HAMAP,PANTHER # (1)!
```

1. CDD, HAMAP, and PANTHER will not be executed


!!! tip

    Analysis names are case-insensitive, and hyphens and underscores are ignored: `CATH-Gene3D`, `cathgene3d`, and `CATH_GENE3D` are all valid.

!!! warning

    Some analyses depend on separately licensed third-party applications. See [Licensed applications](licensed_applications.md).


## Annotate nucleotide sequences

By default, InterProScan expects protein sequences, but supports nucleotide sequences when the `--nucleic` flag is on:

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/contigs.fna \ # (1)!
  --nucleic # (2)!
```

1. Path to a FASTA file of nucleic acid sequences
2. InterProScan will do a six-frame translation before annotatings the translated ORFs

InterProScan uses [Easel](https://github.com/EddyRivasLab/easel) `esl-translate` with its default settings. By default, it reports ORFs longer than 20 amino acids, does not require a specific start codon, and allows any amino acid to begin an ORF. This helps with sequence fragments, genes with introns, and other incomplete coding regions. 

!!! warning

    Because one nucleotide sequence can produce multiple ORFs, nucleotide FASTA files often take much longer to analyse than protein FASTA files of the same size. To improve performance, use stricter ORF filtering before running InterProScan.

## Including GO terms and pathway annotations

InterPro groups related signatures from member databases such as Pfam, PIRSF, and PANTHER into curated InterPro entries. These entries capture protein families, domains, and functional sites, and can carry cross-references such as Gene Ontology (GO) terms and pathway mappings.

When a match is integrated into an InterPro entry, InterProScan can add those annotations to the output.

```bash
nextflow run ebi-pf-team/interproscan6 \
  -r 6.0.0 \
  -profile docker \
  --datadir data \
  --input /path/to/sequences.faa \
  --goterms \ # (1)!
  --pathways # (2)!
```

1. Includes GO terms associated with matched InterPro entries
2. Includes pathway annotations derived from matched InterPro entries

GO terms are curator-assigned to InterPro entries and describe conserved molecular functions, biological processes, or cellular locations. Pathway mappings are inferred from relationships between InterPro entries and reviewed UniProtKB proteins with EC numbers.
