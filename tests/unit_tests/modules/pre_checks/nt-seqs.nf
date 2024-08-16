nextflow_process {

    name "Test Process CHECK_NUCLEIC"
    script "interproscan/modules/pre_checks/main.nf"
    process "CHECK_NUCLEIC"

    stage {
        symlink "interproscan/scripts/pre_checks/check_nucleic_seq.py"
    }

    test("Should pass") {

        when {
            params {
                input="$projectDir/tests/unit_tests/test_inputs/pre_analysis_checks/nt_seqs.fasta"
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
                input="$projectDir/tests/unit_tests/test_inputs/pre_analysis_checks/protein_seqs.fasta"
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
