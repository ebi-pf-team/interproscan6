# Testing

# Unit and integration tests

Unit tests and integration tests are created using [`nf-test`](https://github.com/askimed/nf-test).

The full `nf-test` documentation can be found [here](https://www.nf-test.com/).

All tests are stored in the `tests/` dir.

## Set up

Install using Conda

```
conda install -c bioconda nf-test
```

## Useful arguments

To run tests:
```
nf-test test <path to test file>
```

Useful arguments:
* `--debug` -- Show debugging infos and dump channels
* `--updateSnapshot` -- Use this flag to re-record every snapshot that fails during this test run.
* `--verbose` -- Show Nextflow Ouput -- Very handy when a test assertion is failing
