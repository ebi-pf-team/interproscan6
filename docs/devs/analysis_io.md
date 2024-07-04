# HMMER
## `HMMER_RUNNER` Module
### Generic process to execute HMMER
**Members using:** antifam, ncbifam, gene3d, panther, pfam
- **Input**:
  - Fasta file
  - Runner params 
    - includes: member name, hmm path, switches, member release, postprocessing_params
- **Output**:
  - Path to HMMER `.out` file
  - postprocessing_params

## `HMMER_PARSER` Module
### Generic process to parser HMMER `.out` file
**Members using:** antifam, ncbifam, _funfam_, gene3d, _hamap_, panther, pfam
- **Input**:
  - HMMER `.out` file
  - postprocessing_params
- **Output**:
  - Path to HMMER `.out` file **parsed**


