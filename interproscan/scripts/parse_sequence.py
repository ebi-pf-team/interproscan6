import hashlib
import json
import sys

ILLEGAL_CHARAS = {
    "antifam": "-",
    "cdd": "",
    "coils": "",
    "funfam": "-_.",
    "gene3d": "-_.",
    "hamap": "-_",
    "mobidb": "",
    "ncbifam": "-",
    "panther": "-",
    "pfam": "-",
    "pirsf": "-",
    "prints": "-._",
    "prosite_patterns": "",
    "prosite_profiles": "-._",
    "sfld": "-._",
    "smart": "",
    "superfamily": "-",
    "signalp": "",
    "phobius": "-*._oxuzj",
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


def parse(sequences: dict, applications: str):
    results = {}
    applications = applications.split(",")
    illegal_char_set= {chara for application in applications for chara in ILLEGAL_CHARAS[application]}
    for key, sequence in sequences.items():
        sequence_info = []
        seq = sequence.lower()
        acc = key.split(" ", 1)[0]
        if ">" in sequence:
            raise ValueError(f"{acc} contains illegal character '>'")
        set(seq).intersection(illegal_char_set)
        illegal_matches = set(seq).intersection(illegal_char_set)
        if illegal_matches:
            error_message = f"{acc} contains illegal character(s): \n"
            match_record = {match: [application for application in applications if match in ILLEGAL_CHARAS[application]] for match in illegal_matches}
            for match, tools in match_record.items():
                error_message += f"'{match}' not allowed in {','.join(x for x in tools)} \n"
                #pass
            raise ValueError(error_message)
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
    sequence_parsed = parse(sequences, args[1])
    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
