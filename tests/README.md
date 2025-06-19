# Tests

## Regression tests

* Stored in `./tests/regression_tests`.
* Primary function to validate the structure and content of the final iprscan output files.
* Can also be used to check the content of internal, temporary iprscan files (e.g. the intermediate files created by the `SCAN` subworkflow).

### Requirements

* Python >=3.10

### Usage

```bash
$ python3 ./tests/regression_tests/test_matches.py --help
usage: IPS_match_regression_test [-h] [--expected EXPECTED] [--observed OBSERVED] [--summary] [--format {json,tsv,xml,intermediate}] [--applications APPLICATIONS]

Check presence of matches

options:
  -h, --help            show this help message and exit
  --expected EXPECTED   JSON with expected results (default: tests/data/output/test.faa.json)
  --observed OBSERVED   JSON output file from IPS6 (default: tests/data/output/test.faa.json)
  --summary             Print only the summary message (default: False)
  --format {json,tsv,xml,intermediate}
                        Format of input files. 'json' [default], 'tsv', or 'xml' for final output files or 'intermediate' to compare the temporary working files of InterProScan6. (default: json)
  --applications APPLICATIONS
                        Limit the comparison to a comma-separated list of applicaitons (default: None)
```

## Unit tests

Stored in `./tests/unit_tests`.

### Requirements

* nextflow==24.10.4
* nf-test==0.9.2

### Usage

These test **must** be run from the root of the repository.

To run all tests:
```bash
nt-test test tests/unit_tests
```

To run a specific test point to the specific `nf-test` (`nf.test`) file, for example:
```bash
nf-test test tests/unit_tests/modules/combine/representative_locations/main.nf.test
```
