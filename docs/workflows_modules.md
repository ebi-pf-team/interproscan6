# Workflows, Modules and Executables

This file contains a description of the current workflows, modules and executables in `interproscan-6`. This file aims to introduce and describe the workflows and modules in the order they are implemented in the main (or overall) workflow of `interproscan-6`.

# Main Workflow

The main workflow can be found at `./main.nf`.

Input is a YAML configuration file, which contains the bash parameters from `interproscan-5`.

All parameters (noted by their `params` prefix) are stored in `./nextflow.config`, under `params`, these are:
* `entries`: XML release of InterPro entries
* `goterms`: Gene Ontologoy terms from InterPro release
* `pathways`: ??? 

The values of these parameters can also be configured from their default values (defined in `./nextflow.config`) through the input YAML file.

# Preparing the data for analysis

The first task for `interproscan-6` is to prepare the input data for thedownstream analyses.

## `fasta_channel` Channel

The first step is to initalise the `fasta_channel`.

The channel takes in the path to the input multi-sequence FASTA file, removes duplicate sequences and splits up the file into smaller multi-sequence FASTA files, from here on referred to as "working-FASTA files".

The maximum number of sequences in each of these working-FASTA files is set up be the `batchsize` parameter `params` in `./nextflow.config`.

## `PARSE_SEQUENCE` Module


