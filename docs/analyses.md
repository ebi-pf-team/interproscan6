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

# Analyses

This page lists the analyses available in InterProScan 6 and the names accepted by `--applications` and `--skip-applications`.

!!! tip

    Analysis names are case-insensitive, and hyphens, underscores, and spaces are ignored. For example, `CATH-Gene3D`, `cathgene3d`, and `CATH_GENE3D` are all valid.

## Analysis catalogue

| Analysis | Description | Included by default | Requires separate licensed application |
|----------|-------------|---------------------|----------------------------------------|
| AntiFam | Identifies sequences likely to be spurious or misannotated | Yes | No |
| CATH-Gene3D | Identifies structural domains from the CATH classification | Yes | No |
| CATH-FunFam | Groups protein domains into functional families based on CATH | Yes | No |
| CDD | Detects conserved domains using position-specific scoring matrices from NCBI | Yes | No |
| COILS | Predicts coiled-coil regions based on sequence patterns | Yes | No |
| DeepTMHMM | Predicts transmembrane helices | No | Yes |
| HAMAP | Identifies high-confidence protein families in microbial and organellar proteomes | Yes | No |
| MobiDB-lite | Predicts intrinsically disordered regions | Yes | No |
| NCBIFAM | Matches proteins to curated HMMs from NCBI, including TIGRFAMs | Yes | No |
| PANTHER | Classifies proteins into families and subfamilies with curated GO terms | Yes | No |
| Pfam | Detects protein domains and families using HMMs built from multiple sequence alignments | Yes | No |
| Phobius | Predicts transmembrane topology and signal peptides | No | Yes |
| PIRSF | Classifies proteins into evolutionary families based on full-length sequence similarity | Yes | No |
| PIRSR | Identifies conserved residues using manually curated site rules | Yes | No |
| PRINTS | Detects protein families using groups of conserved motifs | Yes | No |
| PROSITE-patterns | Identifies protein features based on short sequence motifs | Yes | No |
| PROSITE-profiles | Detects protein families and domains using position-specific scoring profiles | Yes | No |
| SFLD | Classifies enzymes by relating sequence features to chemical function | Yes | No |
| SMART | Identifies signaling and extracellular domains | Yes | No |
| SUPERFAMILY | Assigns structural domains using HMMs based on the SCOP superfamily classification | Yes | No |
| SignalP-Euk | Predicts signal peptides in eukaryotic proteins | No | Yes |
| SignalP-Prok | Predicts signal peptides in prokaryotic proteins | No | Yes |
| TMbed | Predicts transmembrane helices, transmembrane strands, and signal peptides | No | No |

## Default selection behaviour

- By default, InterProScan runs all analyses except the deep-learning-based analyses (DeepTMHMM, SignalP-Euk, SignalP-Prok, and TMbed).
- DeepTMHMM, SignalP-Euk, SignalP-Prok, and TMbed are only run when explicitly selected with `--applications`, or when you use `--run-ml`.
- Phobius is not a deep-learning-based analysis, but it still requires separate installation before it becomes available.

## Licensed applications

Some InterProScan analyses rely on separately licensed third-party applications.

This affects the external software component, not your ability to run InterProScan itself or to analyse your own sequences with the bundled analyses.

The licensed analyses are:

- `DeepTMHMM`
- `Phobius`
- `SignalP-Euk`
- `SignalP-Prok`

See [Licensed applications](licensed_applications.md) for installation and activation instructions.

## Changes from InterProScan 5

InterProScan 6 supports most of the analyses available in InterProScan 5.

The following analyses were available in InterProScan 5 but are not
included in InterProScan 6:

| Name      | Reference                                                    | Description                                                              |
|-----------|--------------------------------------------------------------|--------------------------------------------------------------------------|
| SignalP 4 | [Petersen et al., 2011](https://doi.org/10.1038/nmeth.1701)  | Prediction of the presence and location of signal peptide cleavage sites |
| TMHMM     | [Krogh et al., 2001](https://doi.org/10.1006/jmbi.2000.4315) | Prediction of transmembrane helices                                      |

InterProScan 6 also adds the following analyses:

| Name        | Reference                                                         | Description                                                         |
|-------------|-------------------------------------------------------------------|---------------------------------------------------------------------|
| DeepTMHMM   | [Krogh et al., 2001](https://doi.org/10.1101/2022.04.08.487609)   | Prediction of transmembrane helices                                 |
| SignalP 6   | [Teufel et al., 2022](https://doi.org/10.1038/s41587-021-01156-3)   | Prediction of signal peptides and their cleavage sites in all domains of life |
| TMbed       | [Bernhofer & Rost, 2022](https://doi.org/10.1186/s12859-022-04873-x) | Prediction of transmembrane proteins through Language Model embeddings |
