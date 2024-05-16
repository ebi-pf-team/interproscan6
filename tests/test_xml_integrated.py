import pytest
import subprocess
import xml.etree.ElementTree as ET
import os


@pytest.fixture
def get_current_output() -> ET.Element:
    def run_interproscan(input_file, output_file):
        command = f"nextflow run interproscan.nf --input {input_file} --applications antifam --formats xml --disable_precalc --output {output_file}"
        subprocess.run(command, shell=True)
    project_dir = os.path.dirname(os.path.abspath(__file__))
    output_file = "tests/test_outputs/xml_output"
    input_file = "tests/test_outputs/xml_test.fasta"
    run_interproscan(input_file, output_file)
    current_output = project_dir + "/test_outputs/xml_output" + ".xml"
    with open(current_output, 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


@pytest.fixture
def get_expected_output() -> ET.Element:
    expected_output = "tests/test_outputs/xml_output_expected.xml"
    with open(expected_output, 'r') as f:
        tree = ET.parse(f)
    return tree.getroot()


def compare_xml_elements(elem1, elem2):
    if elem1.tag != elem2.tag:
        raise AssertionError(f"Different tags: '{elem1.tag}' != '{elem2.tag}'")

    if elem1.attrib != elem2.attrib:
        raise AssertionError(f"Different attributes: {elem1.attrib} != {elem2.attrib}")

    if (elem1.text or '').strip() != (elem2.text or '').strip():
        raise AssertionError(f"Different texts: '{elem1.text}' != '{elem2.text}'")

    if len(elem1) != len(elem2):
        raise AssertionError(f"Different quantities: {len(elem1)} != {len(elem2)}")

    for child1 in elem1:
        found_match = False
        for child2 in elem2:
            try:
                compare_xml_elements(child1, child2)
                found_match = True
                break
            except AssertionError:
                pass
        if not found_match:
            raise AssertionError(f"Element not found: {child1.tag}")

    return True


def test_xml_output(get_expected_output, get_current_output):
    assert compare_xml_elements(get_expected_output, get_current_output)
