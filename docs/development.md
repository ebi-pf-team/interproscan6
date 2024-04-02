# Development and Observations

This document is for tracking development of `interproscan-6`. This includes implementation of new features, enhancing existing features, and bus fixing.

# Must haves and Must do

### Input

- [ ] Allow user to provide nucleotide sequences

### Documentation

- [x] Add description of configuration file

### Scripts

- [ ] Re-write `ps_scan.pl` in Python in order to remove the perl dependency

# Should haves and Should do

### Performance

- [ ] Test batchsizes - access performance and resource cost

### Inputs

- [ ] Validate input file format - e.g. if user fails to provide sequences in FASTA format, an error message that will be understandable to a none expert will be raised
- [ ] If no sequences are provided, the program should automatically close without initialising the downstream analyses

### Code

- [ ] Parse arguments directly to script, remove reliance on `argsparse`.

This can be done using the `sys` module.

```python
import sys

def my_function(arg1, arg2):
    print(f"Argument 1: {arg1}, Argument 2: {arg2}")

if __name__ == "__main__":
    # sys.argv[0] is the script name itself, so we exclude it.
    args = sys.argv[1:]
    my_function(*args)
```

Then to run this script:

```bash
python my_script.py arg1 arg2
```

### Documentation

- [ ] Add description of workflows, modules and executables (for devs)

# Could haves and Could do

### Configuration

- [ ] Remove `help: false/true` in input config yaml. If no yaml file provided, print help information, else parse input yaml and run workflow.
- [ ] Write separate/additional requirements for dev (typically includes additional linters and doc builders, `requirements-dev.txt`) and running (for users `requirements.txt`)
- [ ] Docs - make note on handling of ambigous amino acids, and accepted symbols
- [ ] Consider building a metamodel
- If `applications` in `input.yaml` is None, run all tools. This will make it easier for users, saving them from having to write all the names of all the tools

# Observations

Food for thoughts - direct actions are necessarily defined but these observations may lead to changes in the future.

## Match lookup

In the lookup we have to filter in the script because when we make the request we don't have the option of passing the applications we want first, the lookup match will return all the applications that have match, only then I can filter (using the `applications` variable).
