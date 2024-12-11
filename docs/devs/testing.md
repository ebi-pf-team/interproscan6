# Testing

Unit tests and integration tests are created using [`nf-test`](https://github.com/askimed/nf-test) [Patel et al., 2024].

> Yash Patel, Chenghao Zhu, Takafumi N Yamaguchi, Yuan Zhe Bugh, Mao Tian, Aaron Holmes, Sorel T Fitz-Gibbon, Paul C Boutros, NFTest: automated testing of Nextflow pipelines, Bioinformatics, Volume 40, Issue 2, February 2024, btae081, https://doi.org/10.1093/bioinformatics/btae081

The full `nf-test` documentation can be found [here](https://www.nf-test.com/).

The in-house python scripts are tested using [`pytest`](https://docs.pytest.org/en/8.2.x/).

## Set up

Install using Conda

```
conda install -c bioconda nf-test
conda install anaconda::pytest
```

## Writing tests

All tests are stored in the `tests/` dir:
* Nextflow unit tests are written in subdirectories
* Python pytests are written in the `tests/` dir

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

To run Nextflow tests:
```
nf-test test <path to test file / or path to dir containing test files>
```

Useful arguments:
* `--debug` -- Show debugging infos and dump channels
* `--updateSnapshot` -- Use this flag to re-record every snapshot that fails during this test run.
* `--verbose` -- Show Nextflow Ouput -- Very handy when a test assertion is failing

To run `pytests` use:
```
python -m pytest tests/unit_tests
```
This method adds the current directory to `sys.path` which is essential for `pytest` to find the 
scripts, owing to the structure of this repository.
