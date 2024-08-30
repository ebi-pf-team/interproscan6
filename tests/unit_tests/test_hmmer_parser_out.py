"""Test the python script scripts/hmmer/test_parser_out.py

These test are intend to be run from the root of the repository using:
pytest -v
"""

import json
import pytest


from interproscan.scripts.hmmer.parser_out import parse


def test_check_output_types(test_input_dir):
    input_path = str(test_input_dir / "pfam_out")
    result = parse(input_path, "pfam")

    sequence_data = result["sp|A2SLW2|1A1D_METPP"]["PF00291"]

    # Ensuring the correct types is extremely important for filtering matches for some members
    assert isinstance(sequence_data["evalue"], float)
    assert isinstance(sequence_data["score"], float)
    assert isinstance(sequence_data["qlen"], int)

    for location in sequence_data["locations"]:
      assert isinstance(location["start"], int)
      assert isinstance(location["end"], int)
      assert isinstance(location["hmmStart"], int)
      assert isinstance(location["hmmEnd"], int)
      assert isinstance(location["hmmLength"], int)
      assert isinstance(location["evalue"], float)
      assert isinstance(location["score"], float)


@pytest.mark.parametrize(
    "input_path, member_db, expected_result",
    [
        ("hmmer3_parser/pfam_out", "pfam", "hmmer3_parser/pfam_parsed.json"),
        ("hmmer3_parser/gene3d_out", "gene3d", "hmmer3_parser/gene3d_parsed.json"),
        ("hmmer3_parser/panther_out", "panther", "hmmer3_parser/panther_parsed.json")
    ]
)
def test_parse_multiple_sequences(test_input_dir, test_output_dir, input_path, member_db, expected_result):
    input_file_path = test_input_dir / input_path
    output_file_path = test_output_dir / expected_result

    result = parse(str(input_file_path), member_db)
    with output_file_path.open('r') as file:
        expected_result = json.load(file)

    assert result == expected_result
