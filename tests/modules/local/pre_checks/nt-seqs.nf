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

    test("NT-seq pre-check should fail") {

        when {
            params {
                input="$projectDir/tests/test_inputs/protein_seqs.fasta"
            }
            process {
                """
                input[0] = file("${params.input}")
                """
            }
        }

        then {
            assert process.failed
            assert process.exitStatus == 1
        }

    }

}
