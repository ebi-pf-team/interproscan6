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

A separate label for each analysis tool runner, including:
* `cdd_runner`
* `coils_runner`
* `hmmer_runner`
* `hmmer_2_runner`
* `mobidb_runner`
* `prosite_pfscan_runner`
* `prosite_pfsearch_runner`
* `signalp_runner`

### Label `analysis_parser`

For all processes in post-processing and parsing the output from the member database analyses, as well as filtering the hits retrieved from parsing the output.

Processes for parsing the output:
* `CDD_PARSER`
* `COILS_PARSER`
* `HMMER_PARSER`
* `HMMER2_PARSER`
* `HMMER_SCAN_PARSER`
* `PFSCAN_PARSER`
* `PFSEARCH_PARSER`
* `MOBIDB_PARSER`
* `SIGNALP_PARSER`

Process for post-processing:
* `ADD_CATH_SUPERFAMILIES`
* `CATH_RESOLVE_HITS`
* `CDD_POSTPROCESS`
* `HAMAP_POST_PROCESSER`
* `PANTHER_POST_PROCESSER`
* `PIRSF_POST_PROCESSER`
* `SFLD_POST_PROCESSER`
* `SUPERFAMILY_POST_PROCESSER`

Processes for filtering matches:
* `FUNFAM_FILTER_MATCHES`
* `GENE3D_FILTER_MATCHES`
* `HAMAP_FILTER_MATCHES`
* `PANTHER_FILTER_MATCHES`
* `PFAM_FILTER_MATCHES`
* `PIRSF_FILTER_MATCHES`
* `SFLD_FILTER_MATCHES`
* `SMART_FILTER_MATCHES`
* `SUPERFAMILY_FILTER_MATCHES`

### Label `get_orfs`

All processes related to the translation of nucleotide sequences into protein sequences.

Processes:
* `GET_ORFS`

### Label `treegrafter_analysis`

The analysis tool `Treegrafter` is used during the post-processing of Panther hits. This process (`PANTHER_POST_PROCESSER`) is given its own 
label to allow us to assign more memory to this specific process, as `Treegrafter` can be memory intensive.

### Label `io`

All processes related to parsing the input files and generating the outputs files.

Input handling processes:
* `PARSE_SEQUENCE`
* `CHECK_NUCLEIC`

Output handling processes:
* `AGGREGATE_RESULTS`
* `WRITE_RESULTS`

Xref processes:


### Label `mls` - Match Lookup Service

All processes related to calling to and parsing data from the Match Lookup Service (MLS).

Processes: 
* `LOOKUP_CHECK`
* `LOOKUP_MATCHES`
* `LOOKUP_NO_MATCHES`

## Label `xref`

All processes related to the XREF subworkflow. These processes retrieve signature and InterPro entry data, adding these data to the hits. If enabled, these processes also retrieve and add GO terms, pathways and (Panther only) PAINT annotation data to the hits.

Processes:

* `ENTRIES`
* `GOTERMS`
* `PAINT_ANNOTATIONS`
* `PATHWAYS`
