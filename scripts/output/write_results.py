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
        for seq_target, info in seq_matches.items():
            sequence_data = info['sequences']
            matches = info["matches"]

            seq_id = seq_target
            md5 = sequence_data[2]
            seq_len = sequence_data[3]
            for match_acc, match in matches.items():
                member_db = match["member_db"]
                sig_acc = match["accession"]
                sig_desc = ""  # info on DB or hmm.out (later step)
                status = "T"
                interpro_acc = "-"  # it will be added on entries PR
                interpro_desc = "-"  # it will be added on entries PR
                for location in match["locations"]:
                    ali_from = location["start"]
                    ali_to = location["end"]

                    tsv_file.write(
                        f"{seq_id}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t{sig_desc}\t{ali_from}\t{ali_to}\t{status}\t{current_date}\t{interpro_acc}\t{interpro_desc}\n")


def json_output(seq_matches: dict, output_path: str):
    json_output = os.path.join(output_path + '.json')
    results = []
    for seq_id, data in seq_matches.items():
        sequence = data['sequences'][1]
        md5 = data['sequences'][2]
        matches = []
        if 'matches' in data and data['matches']:
            for match_key, match_data in data['matches'].items():
                signature = {
                    "accession": match_data['accession'],
                    "description": match_data['name'],
                    "signatureLibraryRelease": {
                        "library": match_data['member_db'].upper(),
                        "version": match_data['version']
                    }
                }

                location = match_data['locations'][0]
                match = {
                    "signature": signature,
                    "locations": [
                        {
                            "start": location['start'],
                            "end": location['end'],
                            "representative": location['representative'],
                            "hmmStart": location['hmmStart'],
                            "hmmEnd": location['hmmEnd'],
                            "hmmLength": location['hmmLength'],
                            "hmmBounds": location['hmmBounds'],
                            "evalue": location['evalue'],
                            "score": location['score'],
                            "envelopeStart": location['envelopeStart'],
                            "envelopeEnd": location['envelopeEnd'],
                            "postProcessed": location['postProcessed']
                        }
                    ],
                    "evalue": match_data['evalue'],
                    "score": match_data['score'],
                    "model-ac": match_key
                }

                matches.append(match)
        result = {
            "sequence": sequence,
            "md5": md5,
            "matches": matches
        }
        results.append(result)

    final_data = {"interproscan-version": "6.0.0", 'results': results}
    with open(json_output, 'w') as json_file:
        json_file.write(json.dumps(final_data, indent=2))


def write_results(sequences_path: str, matches_path: str, output_format: str, output_path: str):
    output_format = output_format.upper()
    seq_matches = {}

    all_sequences = {}
    with open(matches_path, 'r') as match_data:
        all_matches = json.load(match_data)
    with open(sequences_path, 'r') as seq_data:
        for line in seq_data:
            sequence = json.loads(line)
            all_sequences.update(sequence)

    for key in all_sequences:
        if key in all_matches:
            seq_matches[key] = {
                'sequences': all_sequences[key],
                'matches': all_matches[key]
            }
        else:
            seq_matches[key] = {
                'sequences': all_sequences[key],
                'matches': {}
            }

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



# match_data
# "match_data": {
#         "ANF00006": {
#           "accession": "ANF00006",
#           "name": "Spurious_ORF_06",
#           "description": "",
#           "e_value": 1.5e-50,
#           "score": 160.8,
#           "bias": 18.2,
#           "member_db": "AntiFam",
#           "locations": [
#             {
#               "start": 9,
#               "end": 40,
#               "representative": "",
#               "hmmStart": 19,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.5e-09,
#               "score": 28.0,
#               "envelopeStart": 4,
#               "envelopeEnd": 40,
#               "bias": "0.8",
#               "postProcessed": ""
#             },
#             {
#               "start": 52,
#               "end": 101,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.5e-16,
#               "score": 50.4,
#               "envelopeStart": 52,
#               "envelopeEnd": 101,
#               "bias": "0.2",
#               "postProcessed": ""
#             },
#             {
#               "start": 113,
#               "end": 162,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.3e-14,
#               "score": 44.1,
#               "envelopeStart": 113,
#               "envelopeEnd": 162,
#               "bias": "0.2",
#               "postProcessed": ""
#             },
#             {
#               "start": 174,
#               "end": 223,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 3.2e-14,
#               "score": 44.5,
#               "envelopeStart": 174,
#               "envelopeEnd": 223,
#               "bias": "0.1",
#               "postProcessed": ""
#             },
#             {
#               "start": 9,
#               "end": 40,
#               "representative": "",
#               "hmmStart": 19,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.5e-09,
#               "score": 28.0,
#               "envelopeStart": 4,
#               "envelopeEnd": 40,
#               "bias": "0.8",
#               "postProcessed": ""
#             },
#             {
#               "start": 52,
#               "end": 101,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.5e-16,
#               "score": 50.4,
#               "envelopeStart": 52,
#               "envelopeEnd": 101,
#               "bias": "0.2",
#               "postProcessed": ""
#             },
#             {
#               "start": 113,
#               "end": 162,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 4.3e-14,
#               "score": 44.1,
#               "envelopeStart": 113,
#               "envelopeEnd": 162,
#               "bias": "0.2",
#               "postProcessed": ""
#             },
#             {
#               "start": 174,
#               "end": 223,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 50,
#               "hmmLength": 50,
#               "hmmBounds": "",
#               "evalue": 3.2e-14,
#               "score": 44.5,
#               "envelopeStart": 174,
#               "envelopeEnd": 223,
#               "bias": "0.1",
#               "postProcessed": ""
#             }
#           ]
#         },
#         "ANF00057": {
#           "accession": "ANF00057",
#           "name": "Spurious_ORF_57",
#           "description": "",
#           "e_value": 1.8e-105,
#           "score": 338.3,
#           "bias": 29.8,
#           "member_db": "AntiFam",
#           "locations": [
#             {
#               "start": 6,
#               "end": 101,
#               "representative": "",
#               "hmmStart": 16,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.9e-29,
#               "score": 93.0,
#               "envelopeStart": 2,
#               "envelopeEnd": 101,
#               "bias": "4.3",
#               "postProcessed": ""
#             },
#             {
#               "start": 70,
#               "end": 162,
#               "representative": "",
#               "hmmStart": 19,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 3.6e-30,
#               "score": 95.9,
#               "envelopeStart": 69,
#               "envelopeEnd": 162,
#               "bias": "0.8",
#               "postProcessed": ""
#             },
#             {
#               "start": 113,
#               "end": 223,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.8e-36,
#               "score": 115.5,
#               "envelopeStart": 113,
#               "envelopeEnd": 223,
#               "bias": "2.5",
#               "postProcessed": ""
#             },
#             {
#               "start": 174,
#               "end": 243,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 70,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.9e-22,
#               "score": 70.4,
#               "envelopeStart": 174,
#               "envelopeEnd": 243,
#               "bias": "0.9",
#               "postProcessed": ""
#             },
#             {
#               "start": 6,
#               "end": 101,
#               "representative": "",
#               "hmmStart": 16,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.9e-29,
#               "score": 93.0,
#               "envelopeStart": 2,
#               "envelopeEnd": 101,
#               "bias": "4.3",
#               "postProcessed": ""
#             },
#             {
#               "start": 70,
#               "end": 162,
#               "representative": "",
#               "hmmStart": 19,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 3.6e-30,
#               "score": 95.9,
#               "envelopeStart": 69,
#               "envelopeEnd": 162,
#               "bias": "0.8",
#               "postProcessed": ""
#             },
#             {
#               "start": 113,
#               "end": 223,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 111,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.8e-36,
#               "score": 115.5,
#               "envelopeStart": 113,
#               "envelopeEnd": 223,
#               "bias": "2.5",
#               "postProcessed": ""
#             },
#             {
#               "start": 174,
#               "end": 243,
#               "representative": "",
#               "hmmStart": 1,
#               "hmmEnd": 70,
#               "hmmLength": 111,
#               "hmmBounds": "",
#               "evalue": 2.9e-22,
#               "score": 70.4,
#               "envelopeStart": 174,
#               "envelopeEnd": 243,
#               "bias": "0.9",
#               "postProcessed": ""
#             }
#           ]
#         }
#       }
#     },
