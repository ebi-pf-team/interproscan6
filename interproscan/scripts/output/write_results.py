import json
import re
import sys

from format_writer.tsv_output import tsv_pro_output
from format_writer.tsv_output import tsv_output
from format_writer.json_output import json_output
from format_writer.xml_output import xml_output


NT_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=.*$")


def write_results(
    sequences_path: str,
    matches_path: str,
    output_format: list,
    output_path: str,
    version: str,
    nucleic: bool,
):
    seq_matches = {}
    all_sequences = {}
    with open(matches_path, 'r') as match_data:
        all_matches = json.load(match_data)
    with open(sequences_path, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
            all_sequences.update(sequence)

    if not nucleic:
        for key, value in all_sequences.items():
            seq_matches[key] = {
                'sequences': value,
                'matches': all_matches.get(key, {})
            }
    else:
        # map ORF<id> to the original nucleic seq
        for key, value in all_sequences.items():
            # key == orf<id>
            # value eg. == ["orf1 source=Bob coords=143..286 length=48
            # frame=2 desc=", "AAARLRRALARWALKNARPKLRINRRARRWARACCWKNIPPTCCCARL",
            # "f721ded069c8a767dac17d73641c938d", 48]
            result_id = f"{NT_PATTERN.match(value[0]).group(1)}_{key}"
            seq_matches[result_id] = {
                'sequences': value,
                'matches': all_matches.get(key, {})
            }

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path)
    if "TSV-PRO" in output_format:
        tsv_pro_output(seq_matches, output_path)
    if "JSON" in output_format:
        json_output(seq_matches, output_path, version)
    if "XML" in output_format:
        xml_output(seq_matches, output_path, version)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats_str = args[2]
    output_path = args[3]
    version = args[4]
    nucleic = True if args[5] == "true" else False

    formats = set(formats_str.upper().split(','))
    write_results(sequences, matches, formats, output_path, version, nucleic)


if __name__ == "__main__":
    main()
