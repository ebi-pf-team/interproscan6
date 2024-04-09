# Applications

This doc lists the member databases (and associated applications), as well as how the analyses are performed, including parameters, switches, post processing steps and third party requirements.

## Third party tools:

- [`HMMER`](https://academic.oup.com/nar/article-lookup/doi/10.1093/nar/gky448) - Assesses alignment between query sequence and Hidden Markov Models (HMMs)
- [`cath-resolve-hits`](https://doi.org/10.1093/bioinformatics/bty863) - helps to resolve domain matches
- [`TreeGrafter`](https://doi.org/10.1093/bioinformatics/bty625) - a tool for annotating uncharacterised protein sequences, using annotated phylogenetic trees.
- [`rpsblast`](https://www.animalgenome.org/blast/doc/rpsblast.html) - reversed position specific BLAST.
- [`nCoils`](https://doi.org/10.1016/S0076-6879(96)66032-7) - tool for the prediction and analysis of coiled-coil structures
- [`modidb`](https://doi.org/10.1093/nar/gkac1065) - Database and software for intrinsically disordered proteins
- [`phobius`](https://doi.org/10.1016/j.jmb.2004.03.016) - Tool for prediction of transmembrane domains
- [`fingerPRINTScan`](https://doi.org/10.1093/bioinformatics/15.10.799) - Search against FingerPRINTScan with a protein query sequence to identify the closest matching PRINTS sequence motif fingerprints in a protein sequence.
- [`pfsearchV3`](https://doi.org/10.1093/bioinformatics/btt129) - a code acceleration and heuristic to search PROSITE profiles
- [`signalP`](https://www.nature.com/articles/s41587-021-01156-3) - predict signal peptides
- [`TmHMM`](https://doi.org/10.1006/jmbi.2000.4315) - a tool to predict transmembrane domains

# HMMER

This section lists all the tools whose primary analysis is completed using `HMMER`. 

Each member databases requires its own HMMs.

## AntiFam

- `HMMER`

## Gene3D

- `HMMER`
    - _Run inital analyses against HMMs_
- `cath-resolve-hits`
    - _Helps to resolve domain matches_
- `assign_cath_superfamilies`
    - _Takes the output from Cath Resolve Hits and assigns CATH superfamilies based on domain id's_
    - _In house python script_
    Data files:
        - `model_to_family_map.tsv`
        - `discontinuous_regs.pkl.py3`

## FunFam

- `HMMER`
- `cath-resolve-hits`
    - _Helps to resolve domain matches_
- `search.py`
    - _In house python script, parsing `HMMER` and `cath-resolve-hits` output_

## Hamap

- `HMMER`

## NCBIfam

- `HMMER`

## Panther

- `HMMER`

## Pfam

- `HMMER`
- Additional data files:
    - `seed` alignment file `Pfam_A.seed.gz`
    - `clan` file `Pfam-C.gz`
    - `dat` file `Pfam-A.dat.gz`

## PirsF

- `HMMER`

## PirsR

- `HMMER`

## SFLD

- `HMMER`
- Additional data files:
    - site annotation file
    - hierarchy file
- Post-processing achieved using an in house C (`sfld_postprocess.h`, `sfld_postprocess.c`)
    - The executable binary is kept in th `bin/sfld/` directory
    - The C scripts are kept in `scripts/members/sfld/`

## SMART

- `HMMER`

## SUPERFAMILY

- `HMMER`
- Additional data files:
    - `self_hits.tab`
    - `cla.scop`
    - `model.tab`
    - `pdbj95d`

# Member db specific tools

This section discusses member databases where their tool is unqiue to them

## CDD

- `rpsblast`
- data files:
    - CDD library
    - Signature list

## Coils

- `ncoils`

## MobiDB

- `mobidb`

## Phobius

- `phobius`

## PRINTS

- `fingerPRINTScan`

## Prosite

- `pfscanV3`
- `pfsearchV3`
- `pfsearchV3` is coorindated using the in house Python script `run_prosite`
- Also requires the pearl script `ps_scan.pl`, this script should be rewritten in Python to remove the perl dependency

## SignalP

- `signalP`

## TMHMM

- `tmhmm`