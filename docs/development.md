# Development and Observations

This document is for tracking development of `interproscan-6`. This includes implementation of new features, enhancing existing features, and bus fixing.

# Must haves and Must do

### Input

- [ ] Allow user to provide nucleotide sequences

### Documentation

- [X] Add description of configuration file

# Should haves and Should do

- [X] Make applicaiton names case insensitivie for the user
- [ ] Test batchsizes - access performance and resource cost
- [X] Validate input file format - e.g. if user fails to provide sequences in FASTA format, an error message that will be understandable to a none expert will be raised
- [X] If no sequences are provided, the program should automatically close without initialising the downstream analyses
- [X] Parse arguments directly to script, remove reliance on `argsparse`.
- [ ] Add description of workflows, modules and executables (for devs) - IN PROGRESS

# Could haves and Could do

- [X] Remove `help: false/true` in input config yaml. If no yaml file provided, print help information, else parse input yaml and run workflow.
- [ ] Write separate/additional requirements for dev (typically includes additional linters and doc builders, `requirements-dev.txt`) and running (for users `requirements.txt`)
- [ ] Docs - make note on handling of ambigous amino acids, and accepted symbols
- [ ] Consider building a metamodel
- [X] If `applications` in `input.yaml` is None, run all tools. This will make it easier for users, saving them from having to write all the names of all the tools
- [ ] Add singularity support - many clusters do not support docker, so many users may need to create their own Singulatiry images and have to alter the `interproscan6` code base to use Singulatiry. Either hardcode in Singulatiry support along side docker support, or write instructions on how to configure `interproscan6` to use Singulatiry
- [ ] Simplify configuration. Instead of requiring the user to navigate several files and directories to find the write configuration files (e.g. while enabling `SignalP` analyses), confiuration files can be gathered into a single place.
- [ ] Increased acceptance of FASTA file formatting, e.g. protein sequence FASTA downloads from NCBI include spaces between each record, other systems do not.

# Observations

Food for thoughts - direct actions are necessarily defined but these observations may lead to changes in the future.

## Match lookup

In the lookup we have to filter in the script because when we make the request we don't have the option of passing the applications we want first, the lookup match will return all the applications that have match, only then I can filter (using the `applications` variable).
