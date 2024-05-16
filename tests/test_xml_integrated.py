import pytest
import subprocess
import xml.etree.ElementTree as ET
import os


@pytest.fixture
def get_current_output() -> ET.Element:
    def run_interproscan(input_file, output_file):
        command = f"nextflow run interproscan.nf --input {input_file} --applications antifam --formats xml --output {output_file}"
        subprocess.run(command, shell=True)
    project_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = "tests/test_outputs/xml_current_output"
    input_file = "tests/test_outputs/xml_test.fasta"
    run_interproscan(input_file, output_file)
    current_output = project_dir + "/test_outputs/xml_current_output" + ".xml"
    with open(current_output, 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


@pytest.fixture
def get_expected_output() -> ET.Element:
    expected_output = "tests/test_outputs/xml_expected_output.xml"
    with open(expected_output, 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


def xml2dict(element):
    result = {}
    if element.attrib:
        result.update({'attributes': element.attrib})
    for child in element:
        if child.tag not in result:
            result[child.tag] = []
        result[child.tag].append(xml2dict(child))
    for key, value in result.items():
        if isinstance(value, list) and all(isinstance(item, dict) for item in value):
            result[key] = sorted(value, key=lambda x: str(x))
    return result


def compare_dicts(expected, current, ignore_elements):
    for key in expected:
        if key in ignore_elements:
            continue
        if key not in current:
            print(f"Key '{key}' missing in current dict")
            continue
        if isinstance(expected[key], dict):
            compare_dicts(expected[key], current[key], ignore_elements)
        elif isinstance(expected[key], list):
            if len(expected[key]) != len(current[key]):
                print(f"List length mismatch for key '{key}'")
            else:
                for i in range(len(expected[key])):
                    compare_dicts(expected[key][i], current[key][i], ignore_elements)
        else:
            if expected[key] != current[key]:
                print(f"Value mismatch for key '{key}'")
                print(f"  expected: {expected[key]}")
                print(f"  current: {current[key]}")


def test_xml_output(get_expected_output, get_current_output):
    def remove_namespace(element):
        for elem in element.iter():
            if '}' in elem.tag:
                elem.tag = elem.tag.split('}', 1)[1]
            elem.attrib = {k.split('}', 1)[-1].lower(): v for k, v in elem.attrib.items()}

    remove_namespace(get_expected_output)
    remove_namespace(get_current_output)

    expected = xml2dict(get_expected_output)
    current = xml2dict(get_current_output)

    ignore_elements = ['representative', 'hmmer3-location-fragment']
    compare_dicts(expected, current, ignore_elements)
    assert expected == current
