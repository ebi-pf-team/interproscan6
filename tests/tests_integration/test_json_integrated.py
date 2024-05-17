import json
import pytest
import subprocess


def get_current_output(input_path: str, applications: str, disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc}" \
              "--formats json --output current_output"
    subprocess.run(command, shell=True)
    return json.loads(json.dumps("current_output.json", indent=4))


def get_expected_output(expected_output_path: str) -> dict:
    return json.loads(json.dumps(expected_output_path, indent=4))


def compare(expected, current, ignore_elements: list):
    for key in expected:
        if key in ignore_elements:
            continue
        if key not in current:
            print(f"Key '{key}' missing in current dict")
            continue
        if isinstance(expected[key], dict):
            compare(expected[key], current[key], ignore_elements)
        elif isinstance(expected[key], list):
            if len(expected[key]) != len(current[key]):
                print(f"List length mismatch for key '{key}'")
            else:
                for i in range(len(expected[key])):
                    compare(expected[key][i], current[key][i], ignore_elements)
        else:
            if expected[key] != current[key]:
                print(f"Value mismatch for key '{key}'")
                print(f"  expected: {expected[key]}")
                print(f"  current: {current[key]}")


@pytest.mark.parametrize("input_path, output_expected_path, applications, disable_precalc", [])
def test_json_output(input_path: str, output_expected_path: str, applications: str, disable_precalc: bool):
    expected = get_expected_output(output_expected_path)
    current = get_current_output(input_path, applications, disable_precalc)

    ignore_elements = ['representative', 'hmmer3-location-fragment']
    print("Missing elements in current output:")
    compare(expected, current, ignore_elements)  # Check if current has missing elements
    print("Extra elements in current output:")
    compare(current, expected, ignore_elements)  # Check if current has extra elements

    assert expected == current
