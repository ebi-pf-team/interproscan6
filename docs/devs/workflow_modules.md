# Workflows, Modules and Executables

This file contains a description of the current workflows, modules and executables in `interproscan-6`, 
and aims to introduce and describe the workflows and modules in the order they are implemented in the 
`interproscan-6` pipeline.

A visual summary of the workflow is presented in [`interproScan6_workflow.pdf`](https://github.com/ebi-pf-team/interproscan6/blob/main/assets/interproScan6_workflow.pdf).

# Configuration and Main Workflow

The main workflow is located in `./main.nf`.

`InterProScan6` is operated via the commandline, although some utilities are configured via `.config` files.

Parameters (noted by their `params` prefix) are stored in `./conf/xrefs.config`, under `params`, these are:
* `entries`: path to XML release of InterPro entries [REQUIRED]
* `goterms`: path to file containing Gene Ontologoy terms from InterPro release [OPTIONAL]
* `pathways`: path to InterPro release `pathways` file [OPTIONAL]

**The use of `entries` is mandatory in the workflow.**

If the user wants to include the `goterms` or `pathways` information in their outputs,
the paths must be defined in the `./nextflow.config` file. If these paths are not defined, `IPS6` will crash.

Applications and their data files are defined in `./conf/applications.config`.

# Prepare data

The `INIT_PIPELINE` prepares and validates the input data for the downstream analyses. This 
includes checking the necessary data files are present and that the input FASTA file does 
not include illegal characters for he selected member databases.

## Create batches -> `ch_fasta` Channel

The input FASTA file is split into batches, using the `batchSize` defined in `nextflow.config`.

`IPS6` does not remove duplicate sequence, it is up to the user to clean up the input file.

## Prepare the input sequences

If the input sequences are nucleic (i.e. the `--nucleic` flag is used), open reading frames (ORFs)
are predicted using `ESL_TRANSLATE`. Then the ORFs or input protein sequences (i.e. when `--nucleic` is
not used) are parsed into a internal JSON that includes the original sequence, a md5 hash and associated 
sequence meta data.

# Check for pre-calculations

InterPro already contains a series of pre-calculated matches. `interproscan-6` can check the provided
InterPro release to see if any of the input/submitted sequences have already been parsed by InterPro. Where an analysis has been previously performed by InterPro (matches may not always be found), the associated annotation data is retrieved, and the sequence is **not** queried against the models by `interproscan` with the aim to generate _de novo_ annotations.

This operation can be by-passed by using the `--disable_precalc` flag when running `InterProScan-6`.

# Analyse sequences

Calculate matches if there are sequences to be analysed, i.e. if `--disable_precalc` was used or the input 
FASTA file contains sequences not previously analysed by InterPro. This triggers the `SCAN_SEQUENCES` subworkflow
that coordinates running the member database signature recognition methods.

* Configuration:
    * `./conf/applications.config` - define operational parameters to members databases, e.g. binary paths, switches commands, ...
* Input:
    * Sequences to be analysed and the names of the applications to be included in the analysis.
* Executes:
    * The specific modules are dependent on the application but in general this entails:
      * Module `RUNNER` - run the analysis method (e.g. `HMMER`)
      * Module `POST_PROCESS` -- _only some member databases_
      * Module `PARSER` - parse the final output into `Match` objects (defined in `lib`), and can include filtering of the matches.
* Output:
    * A JSON file for each member database per batch

# Aggregate the results

First, gather all the matches for each batch into a single JSON, then aggregate all these JSON files into one JSON.

# Representative domains

Screen the applicable member databases and domains for representative domains.

# Write the output

* JSON
* TSV
* XML
