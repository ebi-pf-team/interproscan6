# HMMER
## RUNNER Module
### `HMMER_RUNNER`
Generic process to execute HMMER
**Members using:** antifam, ncbifam, panther, pfam
- **Input**:
  - Fasta file
  - Runner params 
    - includes: member name, hmm path, switches, member release, postprocessing_params
- **Output**:
  - Path to HMMER `.out` file
  - postprocessing_params

### `HMMER_RUNNER_WITH_ALIGNMENTS`
Process to execute HMMER that requires alignments
**Members using:** sfld, gene3d (funfam requires alignment but has a different process)
- **Input**:
  - Fasta file
  - Runner params 
    - includes: member name, hmm path, switches, member release, postprocessing_params
- **Output**:
  - Path to HMMER `.out` file
  - Path to HMMER `.dtbl` file (used only by sfld. It will be removed after .c script translation to .py)
  - Path to HMMER `.align` file
  - postprocessing_params

### `FUNFAM_HMMER_RUNNER`
FUNFAM has a specific process because Gene3D must be run before FunFam
- **Input**:
  - Fasta file
  - Runner params 
    - includes: member name, hmm path, switches, member release, postprocessing_params, _funfam_cath_superfamilies_
  - applications
- **Output**:
  - Path to HMMER `.out` file
  - postprocessing_params

### `HAMAP_HMMER_RUNNER`
HAMAP has a specific process because it uses HMMER tbl output
- **Input**:
  - Fasta file
  - Runner params 
    - includes: member name, hmm path, switches, member release, postprocessing_params
- **Output**:
  - Path to HMMER `.out` file
  - Path to HMMER `.tbl` file
  - fasta file
  - postprocessing_params


## PARSER Module
### `HMMER_PARSER`
Generic process to parser HMMER `.out` file
**Members using:** antifam, ncbifam, _funfam_, _hamap_, panther, pfam
- **Input**:
  - HMMER `.out` file
  - postprocessing_params
- **Output**:
  - Path to HMMER `.out` file **parsed**

### `HMMER_PARSER_WITH_ALIGNMENTS`
Process to parser HMMER
**Members using:** antifam, ncbifam, _funfam_, _hamap_, panther, pfam
- **Input**:
  - HMMER `.out` file
  - HMMER `.dtbl` file (used only by sfld. It will be removed after .c script translation to .py)
  - postprocessing_params
  - is_sfld (boolean to indicate if parse the .dtbl file and to tell if it needs to retrieve site data)
- **Output**:
  - Path to HMMER `.out` OR `.dtbl` file **parsed**


