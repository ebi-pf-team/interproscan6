import argparse
import json
import os
from datetime import datetime
import xml.etree.ElementTree as ET


def tsv_output(seq_matches: dict, output_path: str, is_pro: bool):
    tsv_output = os.path.join(output_path + '.tsv')
    if is_pro:
        tsv_output = tsv_output + "-pro"
    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        alignment_encoded = ""
        for seq_id, info in seq_matches.items():
            seq_len = info[3]
            md5 = info[2]
            for n_match in range(4, len(info)):
                for domain in info[n_match]["domains"]:
                    acc = domain["accession"].split(".")[0]
                    try:
                        signature_desc = domain["signature_desc"]
                        interpro_acc = domain["interpro_annotations_acc"]
                    except:
                        signature_desc = "-"
                        interpro_acc = "-"
                    ali_from = domain["ali_from"]
                    ali_to = domain["ali_to"]
                    i_evalue = domain["iEvalue"]
                    if is_pro:
                        alignment_encoded = domain["alignment_encoded"]
                    tsv_file.write(f"{seq_id}\t{md5}\t{seq_len}\t{acc}\t{signature_desc}\t{ali_from}\t{ali_to}\t{i_evalue}\t{current_date}\t{interpro_acc}\t{alignment_encoded}\n")


def json_output(seq_matches: dict, output_path: str):
    json_output = os.path.join(output_path + '.json')
    # still need to filter just necessary information for json output
    # concatenated_data = {"interproscan-version": "6.0.0", 'results': data}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(seq_matches, indent=2))


def xml_output(matches: str, output_path: str):
    # A lot of changes! Need to be recreated!
    pass


def write_results(sequences_path: str, matches_path: str, output_format: str, output_path: str):
    output_format = output_format.upper()

    all_sequences = {}
    with open(matches_path, 'r') as match_data:
        all_matches = json.load(match_data)
    with open(sequences_path, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
            all_sequences.update(sequence)
    seq_matches = {key: all_sequences[key] + all_matches[key] for key in all_matches if key != 'null'}

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path, False)
    if "TSV-PRO" in output_format:
        tsv_output(seq_matches, output_path, True)
    if "JSON" in output_format:
        json_output(seq_matches, output_path)
    # if "XML" in output_format:
    #     xml_output(all_matches, all_sequences, output_path)


def main():
    parser = argparse.ArgumentParser(description="Write result file")
    parser.add_argument(
        "-seq", "--sequences", type=str, help="all sequences parsed"
    )
    parser.add_argument(
        "-matches", "--matches", type=str, help="all matches result parsed"
    )
    parser.add_argument("-format", "--format", type=str, help="output format")
    parser.add_argument("-out", "--output_path", type=str, help="output path")

    args = parser.parse_args()
    write_results(args.sequences, args.matches, args.format, args.output_path)


if __name__ == "__main__":
    main()
