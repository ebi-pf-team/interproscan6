import os
import json
import subprocess


def get_current_output(current_output_path: str, input_path: str, applications: str, disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run interproscan.nf --input {input_path} --applications {applications} {disable_precalc} --formats json,tsv-pro --output {current_output_path} --goterms --pathways"
    if os.path.exists(str(current_output_path) + ".json"):
        os.remove(str(current_output_path) + ".json")
    subprocess.run(command, shell=True)
    with open(str(current_output_path) + ".json", 'r') as f:
        return json.load(f)


def get_expected_output(expected_output_path: str) -> dict:
    with open(str(expected_output_path) + ".json", 'r') as f:
        return json.load(f)


def json2dict(obj):
    if isinstance(obj, dict):
        result = {}
        for key, value in obj.items():
            result[key] = json2dict(value)
        return dict(sorted(result.items()))
    elif isinstance(obj, list):
        sorted_list = sorted([json2dict(item) for item in obj], key=lambda x: json.dumps(x))
        if sorted_list and isinstance(sorted_list[0], dict) and 'md5' in sorted_list[0]:
            sorted_list = sorted(sorted_list, key=lambda x: x.get('md5', ''))
        if sorted_list and isinstance(sorted_list[0], dict) and 'model-ac' in sorted_list[0]:
            sorted_list = sorted(sorted_list, key=lambda x: x.get('model-ac', ''))
        return sorted_list
    else:
        return obj


def compare(expected, current, ignore_elements: list):
    for key in expected:
        if key == "md5":
            print(f"md5: {expected[key]}")
        # if key == "accession":
        #     print(f"accession: {expected[key]}")
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
            if str(expected[key]).lower() != str(current[key]).lower():
                print(f"Value mismatch for key '{key}'")
                print(f"  expected: {expected[key]}")
                print(f"  current: {current[key]}")


def test_json_output(input_path, expected_output_path, current_output_path, applications, disable_precalc):
    expected_output = get_expected_output(expected_output_path)
    current_output = get_current_output(current_output_path, input_path, applications, disable_precalc)

    expected = json2dict(expected_output)
    current = json2dict(current_output)

    with open('/Users/lcf/PycharmProjects/interproscan6/tests/tests_integration/temp_expected.json', 'w') as file:
        json.dump(expected, file, indent=2)
    with open('/Users/lcf/PycharmProjects/interproscan6/tests/tests_integration/temp_current.json', 'w') as file:
        json.dump(current, file, indent=2)

    ignore_elements = ['representative']
    compare(expected, current, ignore_elements)
    # compare(current, expected, ignore_elements)

    assert expected == current
