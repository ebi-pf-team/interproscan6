import json
import os
import sys
from datetime import datetime


def tsv_output(seq_matches: dict, output_path: str, is_pro: bool):
    tsv_output = os.path.join(output_path + '.tsv')
    if is_pro:
        tsv_output = tsv_output + "-pro"
    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        alignment_encoded = ""
        for prot_acc, info in seq_matches.items():
            md5 = info[2]
            seq_len = info[3]
            for n_match in range(4, len(info)):
                status = "T"
                member_db = info[n_match]["member_db"]
                sig_acc = info[n_match]["accession"]
                for domain in info[n_match]["domains"]:
                    try:
                        sig_desc = domain["signature_desc"]
                        interpro_acc = domain["interpro_annotations_acc"]
                    except:
                        sig_desc = "-"
                        interpro_acc = "-"
                    ali_from = domain["ali_from"]
                    ali_to = domain["ali_to"]
                    # evalue = domain["iEvalue"]
                    # if is_pro:
                    #     alignment_encoded = domain["alignment_encoded"]
                    tsv_file.write(f"{prot_acc}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t{sig_desc}\t{ali_from}\t{ali_to}\t{status}\t{current_date}\t{interpro_acc}\t{alignment_encoded}\n")


def json_output(seq_matches: dict, output_path: str):
    json_output = os.path.join(output_path + '.json')
    final_data = {"interproscan-version": "6.0.0", 'results': seq_matches}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))


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
    print(json.dumps(seq_matches, indent=4))

    if "TSV" in output_format:
        tsv_output(seq_matches, output_path, False)
    # if "TSV-PRO" in output_format:
    #     tsv_output(seq_matches, output_path, True)
    if "JSON" in output_format:
        json_output(seq_matches, output_path)


def main():
    args = sys.argv[1:]

    sequences = args[0]
    matches = args[1]
    formats = args[2]
    output_path = args[3]

    write_results(sequences, matches, formats, output_path)


if __name__ == "__main__":
    main()
