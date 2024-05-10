# Labels and Groups

Instead of having to list the CPU and and memory requiremetns for each process, we can group processes together using the [`label` directive](https://www.nextflow.io/docs/latest/process.html#label). This allows users to config the computer resource usage depending on their computing architecture and needs.

Make sure each process in `InterProScan6` is given a label:

```groovy
process EXAMPLE_PROC {
    label 'example_label'
}
```

See the Nextflow [docs](https://www.nextflow.io/docs/latest/process.html#label) for naming rules and convensions.

## Labels and groups

This is a list of each of the groups (i.e. processes with the same label) in `InterProScan6`:

### Label `<analysis>_runner`

A separate label for each analysis tool runner, including `hmmer_runner` and `signalp_runner`

Processes:
* `HMMER_RUNNER`
* `SIGNALP_RUNNER`

### Label `analysis_parser`

For all processes in post-processing and parsing the output from the member database analyses.

Processes:
* `HMMER_PARSER`
* `PANTHER_POST_PROCESSER`
* `SFLD_POST_PROCESSER`
* `SIGNALP_PARSER`

### Label `get_orfs`

All processes related to the translation of nucleotide sequences into protein sequences.

Processes:
* `GET_ORFS`

### Label `io`

All processes related to parsing the input files and generating the outputs files.

Processes:
* `PARSE_SEQUENCE`
* `LOOKUP_NO_MATCHES`
* `LOOKUP_MATCHES`
* `LOOKUP_CHECK`
* `ENTRIES`
* `GOTERMS`
* `PATHWAYS`
* `WRITE_RESULTS`
* `AGGREGATE_RESULTS`

