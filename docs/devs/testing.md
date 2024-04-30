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

## Writing tests

### Staging

`InterProScan6` relies heavily on the `$projectDir` parameter, especially when executing scripts 
inside processes. When testing a process that runs a script, ensure to define the symlink to the
respective scripts in the test using the `stage` keyword. For example, to unit test the process
`CHECK_NUCLEIC` in the `subworkflows/pre_checks/main.nf:PRE_CHECKS` subworkflow:

```groovy
nextflow_process {

    name "Test Process CHECK_NUCLEIC"
    script "subworkflows/pre_checks/main.nf"
    process "CHECK_NUCLEIC"

    stage {
        symlink "scripts/pre_checks/check_nucleic_seq.py"
    }

    test("Should pass") {

        when {
            params {
                input="$projectDir/tests/test_inputs/nt_seqs.fasta"
            }
            process {
                """
                input[0] = file("${params.input}")
                """
            }
        }

        then {
            assert process.success
        }

    }
}
```

### Notes

`assert XXX.stdout.contains()` does a straight comparison, it does not see if the output contains a substring, but performs a full string comparison.

## Running tests

To run tests:
```
nf-test test <path to test file>
```

Useful arguments:
* `--debug` -- Show debugging infos and dump channels
* `--updateSnapshot` -- Use this flag to re-record every snapshot that fails during this test run.
* `--verbose` -- Show Nextflow Ouput -- Very handy when a test assertion is failing
