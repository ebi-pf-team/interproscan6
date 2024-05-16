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


def compare_xml_elements(elem1, elem2, path=""):
    def remove_namespace(element):
        for elem in element.iter():
            if '}' in elem.tag:
                elem.tag = elem.tag.split('}', 1)[1]
            elem.attrib = {k.split('}', 1)[-1].lower(): v for k, v in elem.attrib.items()}

    remove_namespace(elem1)
    remove_namespace(elem2)

    if elem1.tag.lower() != elem2.tag.lower():
        raise AssertionError(f"Different tags at {path}: '{elem1.tag}' != '{elem2.tag}'")

    if sorted(elem1.attrib.items()) != sorted(elem2.attrib.items()):
        raise AssertionError(f"Different attributes at {path}: {elem1.attrib} != {elem2.attrib}")

    if (elem1.text or '').strip().lower() != (elem2.text or '').strip().lower():
        raise AssertionError(f"Different texts at {path}: '{elem1.text}' != '{elem2.text}'")

    if len(elem1) != len(elem2):
        raise AssertionError(f"Different number of children at {path}: {len(elem1)} != {len(elem2)}")

    if not all(isinstance(child, ET.Element) for child in elem1):
        return True

    elem1_children = sorted(elem1, key=lambda x: x.tag.lower() + (x.get('md5') or '').lower() + (x.get('ac') or '').lower())
    elem2_children = sorted(elem2, key=lambda x: x.tag.lower() + (x.get('md5') or '').lower() + (x.get('ac') or '').lower())

    for child1, child2 in zip(elem1_children, elem2_children):
        compare_xml_elements(child1, child2, path=f"{path}/{child1.tag}")

    return True


def test_xml_output(get_expected_output, get_current_output):
    assert compare_xml_elements(get_expected_output, get_current_output)
