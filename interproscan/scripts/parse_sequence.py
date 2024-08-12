import hashlib
import json
import sys

ILLEGAL_CHARAS = {
    "antifam": ["-"],
    "cdd": [],
    "coils": [],
    "funfam": ["-", "*", "_", "."],
    "gene3d": ["-", "*", "_", "."],
    "hamap": ["-", "*", "_", "b", "o", "x", "u", "z", "j"],
    "mobidb": [],
    "ncbifam": ["-"],
    "panther": ["-", "*"],
    "pfam": ["-", "*"],
    "pirsf": ["-"],
    "prints": ["-", ".", "_"],
    "prosite_patterns": [],
    "prosite_profiles": ["-", ".", "_", "*"],
    "sfld": ["-", ".", "_"],
    "smart": [],
    "superfamily": ["-"],
    "signalp": [],
    "phobius": ["-", "*", ".", "_", "o", "x", "u", "z", "j"],
}


def is_fasta_check(fasta_file: str):
    with open(fasta_file, "r") as f:
        line = f.readline().strip()
        if not line.startswith(">"):
            raise ValueError(f"{fasta_file} is not in FASTA format")


def get_sequences(fasta_file: str) -> dict:
    sequences = {}
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                seq_key = line[1:]
                sequences[seq_key] = ""
            else:
                sequences[seq_key] += line
    return sequences


def check_sequence(sequences: dict, applications: str):
    illegal_char_list = set()
    applications = applications.split(",")
    for application in applications:
        for chara in ILLEGAL_CHARAS[application]:
            illegal_char_list.add(chara)

    for key, sequence in sequences.items():
        if ">" in sequence:
            raise ValueError(f"{key} contains illegal character '>'")
        for i in illegal_char_list:
            if i in sequence.lower():
                app_list = [
                    application
                    for application in applications
                    if i in ILLEGAL_CHARAS[application]
                ]
                raise ValueError(
                    f"{key} contains illegal character '{i}' "
                    f"which cannot be used with {','.join(x for x in app_list)}"
                )

    return sequences


def parse(sequences: dict):
    results = {}
    for key, sequence in sequences.items():
        sequence_info = []
        acc = key.split(" ", 1)[0]
        sequence_info.append(key)
        sequence_info.append(sequence)
        sequence_info.append(hashlib.md5(sequence.encode()).hexdigest())
        sequence_info.append(len(sequence))
        results[acc] = sequence_info
    return results


def main():
    args = sys.argv[1:]
    is_fasta_check(args[0])
    sequences = get_sequences(args[0])
    sequences = check_sequence(sequences, args[1])
    sequence_parsed = parse(sequences)
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
