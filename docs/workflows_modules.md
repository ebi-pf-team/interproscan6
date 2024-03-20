# Workflows, Modules and Executables

This file contains a description of the current workflows, modules and executables in `interproscan-6`, and aims to introduce and describe the workflows and modules in the order they are implemented in the `interproscan-6` pipeline.

# Configuration and Main Workflow

The main workflow can be found in `./main.nf`.

Input is a YAML configuration file, which contains the bash parameters from `interproscan-5`.

All parameters (noted by their `params` prefix) are stored in `./nextflow.config`, under `params`, these are:
* `entries`: path to XML release of InterPro entries
* `goterms`: path to file containing Gene Ontologoy terms from InterPro release
* `pathways`: path to InterPro release `pathways` file

The use of `entries` is mandatory in the workflow, but if the user wants to have the `goterms` or `pathways` informations, it is necessary indicate it in the input YAML by setting the respective values to "true".

# Prepare data

The first task for `interproscan-6` is to prepare the input data for the downstream analyses.

## `fasta_channel` Channel

The first step is to initalise the `fasta_channel`.

* Input:
    * Path to multi-sequence FASTA file
* Executes:
    * Removes duplicates
    * Splits input file into multiple, smaller FASTA files
* Output:
    * Paths to FASTA files

The maximum number of sequences in each of these working-FASTA files is set up be the `batchsize` parameter `params` in `./nextflow.config`.

## `PARSE_SEQUENCE` Module

* Input:
    * Path to FASTA files (from the `fasta_channel`)
* Executes:
    * Python script `scripts/parse_sequence.py`
        * Hashes (MD5) sequences (including their metadata: ID, desc, etc)
* Output:
    * `JSON` file of hashed sequences

# Check for pre-calculations

InterPro already contains a series of pre-calculated matches. `interproscan-6` can check the provided InterPro release to see if any of the input/submitted sequences have already been parsed by InterPro. Where an analysis has been previously performed by InterPro (matches may not always be found), the associated annotation data is retrieved, and the sequence is **not** queried against the models by `interproscan` with the aim generate _de novo_ annotations.

This operation can be by-passed by setting `disable_precalc` to `true` in the configuration YAML file.

## `sequence_precalc` Subworkflow

The subworkflow is defined in `subworkflows/sequence_precalc/main.nf`.

**Configuration:**
The subworkflow is configured using the `subworkflows/sequence_precalc/lookup.config` file, which is pre-populated with the URL to the InterPro match-lookup service and its slugs.

**Input**:
* Path to `JSON` file containing the hashed sequences
* Applications (member databases) the user wants

**Modules:**
The subworkflows incorporates three modules (in order):
1. `LOOKUP_CHECK`
    * Input:
        * Hashed sequences
    * Executes:
        * Python script `scripts/lookup/lookup_check.py`
            * Look up to see if there are any pre-calculated matches from any member database in InterPro
    * Output:
        * `dict` of seqs with matches in InterPro, seqs without matches in InterPro, and all seq data
2. `LOOKUP_MATCHES`
    * Input:
        * `dict` from `LOOKUP_CHECK`
        * List of applications to use in analysis
    * Executes:
        * Python script `scripts/lookup/lookup_matches.py`
            * Retrieves pre-calculated match data for only the specified application
    * Output: 
        * `JSON` dump of parsed matches
3. `LOOKUP_NO_MATCHES`
    * Input:
        * `dict` from `LOOKUP_CHECK`
    * Executes:
        * Python script `scripts/lookup/lookup_no_matches.py`
            * Writes out FASTA seqs of hashed seqs not pre-calculated by InterPro
    * Output:
        * FASTA file of sequences to be analysed by InterProScan (`no_match_lookup_fasta.fasta`)

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

Calculate matches if there are sequences to be analysed, i.e. if `sequence_precalc` was disabled or the input FASTA file contains sequences not previously analysed by InterPro.

## Combine `applications_channel` and (`fasta_channel` OR `sequences_to_analyse`)

**Input**: All sequences that were identified as having not been previously analysed by InterPro (or all them, in case of `disable_precalc`) and combines these with an `applications_channel` to get a cartesian product of all applications we will analyse and the subsets of fasta files. 

This way we parallelize the workflow as `number_of_applications x number_of_splitted_fasta_files` flows. You can better visualize what happens in the example below:

![image](https://github.com/ebi-pf-team/interproscan6/assets/17861151/7310f97d-cec3-4d63-8c13-a399a5fb9ef4)


### `SEQUENCE_ANALYSIS` Subworkflow

If `input_yaml.disable_precalc` is true, and/or there are sequences to analyse following checking for precalculated matches, the module `SEQUENCE_ANALYSIS` is used to coordinate checking for matches against the user specified applications (i.e. member databases).

* Configuration:
    * `subworkflows/sequence_analysis/members.config` - define operational parameters to members databases, e.g. binary paths, switches commands, ...
* Input:
    * Sequences to be analysed and the names of the applications to be included in the analysis.
    * `TSV_PRO`: true/false to generate a .tsv-pro file (adding the cigar alignment in the analysis)
* Executes:
    * Module `RUNNER`
    * Module `PARSER`
* Output:
    * Parsed output from sequence analysis tools used for each members

-----------------------------------------------------------------------------------------------------------------------------
PS: This block is being refactored to be more generic and to support different analysis tools (some members don't use hmmer).

#### `HMMER_RUNNER` Module

* Input:
    1. FASTA sequences to be analysed
    2. Path to HMM profiles
    3. `switches`: operational arguments, e.g. the number cpus to use and thresholds for matches
* Executes:
    * `HMMer`
* Output: 
    1. Path to `HMMer` `.out` file
    2. Path to `HMMer` `.dtbl` file

#### `HMMER_PARSER` Module

* Input:
    1. The output from `HMMER_RUNNER` (path to the `HMMer` output files)
    2. `tsv_pro` - ????
* Executes:
    * If `tsv_pro` is true: `hmmer_parser_out`
    * If `tsv_pro` is false: `hmmer_parser_domtbl`
    * Both scripts parse the output from `HMMer` into  `JSON` file
* Output:
    * `hmmer_parsed_<output>.json`
 
-----------------------------------------------------------------------------------------------------------------------------

# Build Output

Compile and write the outputs from the precalculated analyses from InterPro releases with the matches calculated during the sequence analysis stage.

## `UNION_RESULTS` Module

* Input:
    * Pre-calculated results
    * Output from analysis
* Executes:
    * `scripts/post_proc/union_results.py`
        * Combines results into a `JSON` file
* Output:
    * `JSON` file

## `XREFS` Module

* Input:
    * Path to all matches resulted from the previous steps
    * Path to InterPro entries file
    * Path to goterms file (if `true` in input.yaml)
    * Path to pathways file (if `true` in input.yaml)
* Executes:
    * `scripts/xrefs.py`
* Output:
    * Path to file containing the matches with the xref results

## `WRITE_RESULTS` Module

* Input:
    * All sequences
    * All matches
    * Output file formats - defined by user at the configuration stage (e.g. json, xml, tsv)
    * Path for output files
* Executes:
    * Collect input sequences
    * `scripts/write_output.py`
* Output:
    * One file per output file extension
