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


def compare_xml_elements(elem1, elem2, path=""):
    def remove_namespace(element):
        for elem in element.iter():
            if '}' in elem.tag:
                elem.tag = elem.tag.split('}', 1)[1]
            elem.attrib = {k.split('}', 1)[-1]: v for k, v in elem.attrib.items()}

    def get_protein_name(element):
        xref = element.find("xref")
        if xref is not None:
            return xref.get("name", "")
        return ""

    remove_namespace(elem1)
    remove_namespace(elem2)
    elem1_attrib = {k.lower(): v for k, v in elem1.attrib.items()}
    elem2_attrib = {k.lower(): v for k, v in elem2.attrib.items()}
    elem1_text = (elem1.text or '').strip().lower()
    elem2_text = (elem2.text or '').strip().lower()

    current_path = f"{path}/{elem1.tag}"

    if elem1.tag != elem2.tag:
        print(f"Different tags at {current_path}: '{elem1.tag}' != '{elem2.tag}'")
        return False

    if elem1_attrib != elem2_attrib:
        print(f"Different attributes at {current_path}: {elem1_attrib} != {elem2_attrib}")
        return False

    if elem1_text != elem2_text:
        print(f"Different texts at {current_path}: '{elem1_text}' != '{elem2_text}'")
        return False

    if len(elem1) != len(elem2):
        print(f"Different number of children at {current_path}: {len(elem1)} != {len(elem2)}")
        return False

    for child1 in elem1:
        found_match = False
        for child2 in elem2:
            if compare_xml_elements(child1, child2, path=current_path):
                found_match = True
                break
        if not found_match:
            protein_name = get_protein_name(elem1)
            print(f"Element not found at {current_path}: {child1.tag}, protein name: {protein_name}")

    return True


def test_xml_output(get_expected_output, get_current_output):
    assert compare_xml_elements(get_expected_output, get_current_output)
