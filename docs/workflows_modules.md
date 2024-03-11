# Workflows, Modules and Executables

This file contains a description of the current workflows, modules and executables in `interproscan-6`. This file aims to introduce and describe the workflows and modules in the order they are implemented in the `interproscan-6` pipeline.

# Main Workflow

The main workflow can be found in `./main.nf`.

Input is a YAML configuration file, which contains the bash parameters from `interproscan-5`.

All parameters (noted by their `params` prefix) are stored in `./nextflow.config`, under `params`, these are:
* `entries`: XML release of InterPro entries
* `goterms`: Gene Ontologoy terms from InterPro release
* `pathways`: ??? 

The values of these parameters can also be configured from their default values (defined in `./nextflow.config`) through the input YAML file.

# Prepare data

The first task for `interproscan-6` is to prepare the input data for thedownstream analyses.

## `fasta_channel` Channel

The first step is to initalise the `fasta_channel`.

The channel takes in the path to the input multi-sequence FASTA file, removes duplicate sequences and splits up the file into smaller multi-sequence FASTA files, from here on referred to as "working-FASTA files".

The maximum number of sequences in each of these working-FASTA files is set up be the `batchsize` parameter `params` in `./nextflow.config`.

## `PARSE_SEQUENCE` Module

* Input: path to FASTA file of sequences to be analysed
* Executes: python script `scripts/parse_sequence.py`
    * Loads sequences into memory
    * hashes (MD5) sequences (including their metadata: ID, desc, etc)
* Output: JSON file of hashed sequences

# Check for pre-calculations

InterPro already contains a series of pre-calculated matches. `interproscan-6` can check the provided InterPro release to see if any of the input/submitted sequences have already been parsed by InterPro. Where an existing match is found, the associated annotation data is retrieved and the sequence is not queried against the models by `interproscan` to generate _de novo_ annotations.

This operation can be by-passed by `disable_precalc` in the configuration YAML file to "true".

## `sequence_precalc` Subworkflow

The subworkflow is defined in `subworkflows/sequence_precalc/main.nf`.

The subworkflow is configured using the `subworkdlows/sequence_precalc/lookup.config` file, which is pre-populated with the URL to the InterPro match-lookup service, and paths to directories.... ????

**Input**:
* Path to JSON file containing the hashed sequences

**Modules:**

The subworkflows incorporates three modules (in order):
1. `LOOKUP_CHECK`
    * Input: hashed sequences
    * Executes: Python script `scripts/lookup/lookup_check.py`
        * Look up if pre-calculated matches from any member database in InterPro
    * Output: `dict` of seqs with matches in InterPro, seqs without matches in InterPro, and all seq data
2. `LOOKUP_MATCHES`
    * Input:
        * `dict` from `LOOKUP_CHECK`
        * List of applications to use in analysis
    * Executes: Python script `scripts/lookup/lookup_matches.py`
        * Retrieves pre-calculated match data for only the specified application
    * Output: 
3. `LOOKUP_NO_MATCHES`
    * Input: `dict` from `LOOKUP_CHECK`
    * Executes: scripts/lookup/lookup_no_matches.py
        * Writes out FASTA seqs of hashed seqs where 
    * Output: FASTA file of sequences to be analysed by InterProScan (`no_match_lookup_fasta.fasta`)

**Note:**

Sometimes a sequence has been analysed during the InterPro release process and no matches or sites were found. This is still counted as pre-calculated matches/sites.

The results for an MD5 that is not in InterPro and an MD5 hash that is in InterPro but no matches were found in the last release are identical.

For example, `https://www.ebi.ac.uk/interpro/match-lookup/matches/?md5=SOMEMD5WEDONTHAVEINOURDATABASE` returns:

```xml
<kvSequenceEntryXML>
<matches/>
</kvSequenceEntryXML>
```

`LOOKUP_CHECK` is used to differentiate between these two cases.

**Output:**

* Sequences to be analysed because pre-calculated matches are not available
* Retrieved pre-calculated match and site data

# Analyse sequences

Calculate matches if there are sequences to be analysed (if `sequence_precalc` was disabled or the input FASTA file contains sequences not previously analysed by InterPro).

## `applications_channel` Channel

## `SEQUENCE_ANALYSIS` Module

# Build Output

...

## `UNION_RESULTS` Module

## `XREFS` Module

## `WRITE_RESULTS` Module


