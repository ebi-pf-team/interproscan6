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

# Migrating from InterProScan 5

InterProScan 5 was introduced in 2014 ([Jones et al., 2014](https://doi.org/10.1093/bioinformatics/btu031)). Its output formats reflected the structure and constraints of that release and became correspondingly complex over time. InterProScan 6 preserves the same core biological concepts, but it does not reproduce the InterProScan 5 output formats exactly. Treat it as a migration rather than a drop-in replacement, especially if you have downstream parsers, reports, or workflows that depend on the exact InterProScan 5 layout.

## Command line options

The tables below are a migration guide, not a promise that old command lines will run unchanged. Some option names map cleanly; others were split, renamed, or dropped because InterProScan 6 has a different execution model.

### Input control

| InterProScan 5           | InterProScan 6      |
|--------------------------|---------------------|
| `-i, --input`            | `--input`           |
| `-t, --seqtype`          | `--nucleic`         |

!!! info

    In InterProScan 5, `-t, --seqtype` took a value:
    `p` for proteins (the default) or `n` for nucleotides.
    In InterProScan 6, that choice became simpler:
    add `--nucleic` for nucleotide input, or omit it for protein input.


### Execution and resource control

| InterProScan 5                   | InterProScan 6        |
|----------------------------------|-----------------------|
| `-appl, -applications`           | `--applications`      |
| `-cpu, --cpu`                    | `--max-workers`       |
| `-dp, --disable-precalc`         | `--no-matches-api`    |
| `-exclappl, --excl-applications` | `--skip-applications` |

!!! tip

    Do not translate InterProScan 5 `--cpu` mechanically.
    InterProScan 6 also has a `--cpus` option, but it means something different.
    `--cpus` controls the number of CPUs assigned to each task, whereas
    `--max-workers` controls the maximum number of tasks that can run in parallel
    when InterProScan 6 is executed locally.
    `--max-workers` has no effect when tasks are submitted to a compute environment
    such as an HPC scheduler or cloud provider.

### Output control

| InterProScan 5           | InterProScan 6      |
|--------------------------|---------------------|
| `-b, -output-file-base`  | `--outprefix`       |
| `-d, --output-dir`       | `--outdir`          |
| `-f, --formats`          | `--formats`         |
| `-goterms, --goterms `   | `--goterms`         |
| `-o, --outfile`          | N/A                 |
| `-pa, --pathways`        | `--pathways`        |
| `-T, --tempdir`          | `-w, -work-dir`     |

!!! warning

    InterProScan 5 output naming had multiple overlapping controls.
    InterProScan 6 separates them more cleanly, so old assumptions about exact
    filenames usually need to be revisited.

    In InterProScan 5:

    - `-o, --outfile` specifies the *exact* path of the output file, so it can only
      be used when a single output format is selected with `-f, --formats`.
    - `-b, -output-file-base` and `-d, --output-dir` are mutually exclusive.
      You can either specify the base path of the output file
      (including directories, for example `path/to/output/file`), in which case
      the output extension(s) are added automatically, or specify the output
      directory, in which case the input FASTA filename is used and the
      extension(s) are added automatically.

    In InterProScan 6:

    - `--outdir` controls the *directory* where output files are created.
    - `--outprefix` controls the *base name* of the output files, and the
      extension is added automatically.
    - There is no direct replacement for "write exactly this one filename".
      Build the expected filename from `--outdir`, `--outprefix`, and the selected format.
    - For instance, with `--outdir results --outprefix my_proteome --formats gff3,json` the following files are created:
        - `results/my_proteome.gff3`
        - `results/my_proteome.json`

!!! info

    `-w, --work-dir` is a **Nextflow** option. By default, Nextflow creates a
    `work` directory in the current working directory. Unlike InterProScan 5,
    which cleaned up transient data after the run completed, Nextflow keeps the
    working data for resume and debugging unless you remove it explicitly.


## Analyses

Most InterProScan 5 analyses are still available in InterProScan 6, but not all of them. If your workflow depended on a specific tool rather than just the general annotation outcome, check that dependency explicitly.

The following analyses were available in InterProScan 5 but are not included in InterProScan 6:

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

!!! info

    SignalP supports two analysis modes: one for prokaryotic sequences and one for eukaryotic sequences.
    If you previously treated SignalP as a single fixed analysis, update that assumption in your pipeline configuration.

## Output files

The good news is that the high-level data model is still recognisable: `protein -> matches -> locations`. The less convenient part is that the serialised formats are different enough that old XPath, JSONPath, field-name, and exact-file comparisons should be reviewed rather than trusted.

The main differences are:

- InterProScan 6 adds a top-level `interpro-version` field/attribute alongside `interproscan-version`.
- XML no longer uses analysis-specific element names. It switches to a generic `match/location/fragment/site` structure.
- JSON adds a `source` field on each match. In the current output it is the same as `signature.signatureLibraryRelease.library`, but it is intended to carry the origin of a match in future model-based outputs, while `signature.signatureLibraryRelease.library` continues to identify the target member database.
- Some library names are normalised in InterProScan 6, so literal string checks may need updating:

| InterProScan 5 | InterProScan 6 |
|---|---|
| `CDD` | `CDD` |
| `COILS` | `COILS` |
| `FUNFAM` | `CATH-FunFam` |
| `GENE3D` | `CATH-Gene3D` |
| `HAMAP` | `HAMAP` |
| `MOBIDB_LITE` | `MobiDB-lite` |
| `NCBIFAM` | `NCBIFAM` |
| `PANTHER` | `PANTHER` |
| `PFAM` | `Pfam` |
| `PIRSF` | `PIRSF` |
| `PIRSR` | `PIRSR` |
| `PRINTS` | `PRINTS` |
| `PROSITE_PATTERNS` | `PROSITE patterns` |
| `PROSITE_PROFILES` | `PROSITE profiles` |
| `SMART` | `SMART` |
| `SUPERFAMILY` | `SUPERFAMILY` |

- Some per-location fields are now exposed on the generic `location` object, for example `alignment`, `cigarAlignment`, `level` and `sites`.

For example, InterProScan 5 XML encoded the analysis in the element name:

```xml
<hmmer3-match evalue="5.5E-176" score="598.2">
  <signature ac="PF00183" name="HSP90" type="FAMILY">
    <signature-library-release library="PFAM" version="38.1"/>
  </signature>
</hmmer3-match>
```

InterProScan 6 XML uses a generic `match` element and records the origin in `source`:

```xml
<match source="Pfam" evalue="5.5E-176" score="598.2">
  <signature ac="PF00183" name="HSP90" type="Family">
    <signature-library-release library="Pfam" version="38.1"/>
  </signature>
</match>
```

The same simplification shows up in JSON:

```json
{
  "source": "PROSITE patterns",
  "signature": {
    "accession": "PS00298",
    "signatureLibraryRelease": {
      "library": "PROSITE patterns"
    }
  }
}
```

When migrating downstream pipelines:

1. In XML, replace selectors based on analysis-specific tags such as `hmmer3-match` with selectors over generic `match/location/fragment/site` elements.
2. Use `signature.signatureLibraryRelease.library` to identify the member database, and treat `source` as the match-origin field.
3. Update hard-coded library names to the InterProScan 6 labels.
4. Change tests to compare InterProScan 5 and 6 outputs semantically rather than expecting byte-for-byte equality.
