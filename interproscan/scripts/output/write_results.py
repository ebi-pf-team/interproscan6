import json
import os
import re
import sys

try:   # Nextflow uses the relative imports
    from format_writer.tsv_output import tsv_pro_output
    from format_writer.tsv_output import tsv_output
    from format_writer.json_output import build_json_output_nucleic, build_json_output_protein
    from format_writer.xml_output import build_xml_output_nucleic, build_xml_output_protein
except ModuleNotFoundError:  # but pytest needs the imports from interproscan
    from interproscan.scripts.output.format_writer.tsv_output import tsv_pro_output
    from interproscan.scripts.output.format_writer.tsv_output import tsv_output
    from interproscan.scripts.output.format_writer.json_output import build_json_output_nucleic, build_json_output_protein
    from interproscan.scripts.output.format_writer.xml_output import build_xml_output_nucleic, build_xml_output_protein


NT_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=.*$")


def write_results(
    sequences_path: str,
    matches_path: str,
    output_format: list,
    out_file_name: str,
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
            result_id = f"{NT_PATTERN.match(value['seq_id']).group(1)}_{key}"
            seq_matches[result_id] = {
                'sequences': value,
                'matches': all_matches.get(key, {})
            }

    if "TSV" in output_format:
        tsv_output(seq_matches, os.path.join(out_file_name + '.ips6.tsv'))
    if "TSV-PRO" in output_format:
        tsv_pro_output(seq_matches, os.path.join(out_file_name + '.ips6.tsv-pro.tsv'))
    if "JSON" in output_format:
        if nucleic:
            build_json_output_nucleic(seq_matches, os.path.join(out_file_name + '.ips6.json'), version)
        else:
            build_json_output_protein(seq_matches, os.path.join(out_file_name + '.ips6.json'), version)
    if "XML" in output_format:
        if nucleic:
            build_xml_output_nucleic(seq_matches, os.path.join(out_file_name + '.ips6.xml'), version)
        else:
            build_xml_output_protein(seq_matches, os.path.join(out_file_name + '.ips6.xml'), version)


def main():
    """CL input:
    0. str repr of path to sequences json output file
    1. str repr of path to matches json output file
    2. str repr of output formats separated by commas
    3. str repr of output file name
    4. str repr of InterProScan version
    5. str repr of nucleic acid sequences flag"""

    args = sys.argv[1:]
    sequences = args[0]
    matches = args[1]
    formats_str = args[2]
    out_file_name = args[3]
    version = args[4]
    nucleic = True if args[5] == "true" else False

    formats = set(formats_str.upper().split(','))
    write_results(sequences, matches, formats, out_file_name, version, nucleic)


if __name__ == "__main__":
    main()
