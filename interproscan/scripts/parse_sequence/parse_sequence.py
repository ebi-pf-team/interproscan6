import hashlib
import json
import re
import sys


NT_SEQ_ID_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=(\d+\.+\d+)\s+.+frame=\d+\s+desc=.*$")
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
    "pirsr": "-._*",
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


def get_sequences(fasta_file: str) -> dict[str, str]:
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


def get_nucleic_seqs(fasta_file: str) -> dict[str, dict]:
    sequences = {}
    with open(fasta_file, "r") as f:
        lines = f.readlines()
        for line in lines:
            line = line.strip()
            if line.startswith(">"):
                # split on the first white space to match GETORF handling of seqIDs
                seq_key = line[1:].split(maxsplit=1)[0]
                sequences[seq_key] = {'sequence': "", "name": line[1:]}
            else:
                sequences[seq_key]['sequence'] += line
    return sequences

def parse(sequences: dict, applications: str) -> dict[str, list]:
    results = {}
    applications = applications.split(",")
    illegal_char_set = {chara for application in applications for chara in
                        ILLEGAL_CHARAS[application]}
    for key, sequence in sequences.items():
        sequence_info = []
        seq = sequence.lower()
        acc = key.split(" ", 1)[0]
        if ">" in sequence:
            raise ValueError(f"{acc} contains illegal character '>'")
        illegal_matches = set(seq).intersection(illegal_char_set)
        if illegal_matches:
            error_message = f"{acc} contains illegal character(s): \n"
            match_record = {match: [application for application in applications if match in ILLEGAL_CHARAS[application]] for match in illegal_matches}
            for match, tools in match_record.items():
                error_message += f"'{match}' not allowed in {','.join(x for x in tools)}\n"
            print(error_message, file=sys.stderr)
            sys.exit(1)
        sequence_info.append(key)
        sequence_info.append(sequence)
        sequence_info.append(hashlib.md5(sequence.encode()).hexdigest())
        sequence_info.append(len(sequence))
        results[acc] = sequence_info
    return results


def parse_nucleic(sequences: dict, nt_seqs: dict) -> dict[str, list]:
    results = {}
    for key, sequence in sequences.items():
        sequence_info = []
        acc = key.split(" ", 1)[0]
        sequence_info.append(key)
        sequence_info.append(sequence)
        sequence_info.append(hashlib.md5(sequence.encode()).hexdigest())
        sequence_info.append(len(sequence))
        nt_seq_id = NT_SEQ_ID_PATTERN.match(key).group(1)
        sequence_info.append(nt_seqs[nt_seq_id]['name'])
        sequence_info.append(hashlib.md5(nt_seqs[nt_seq_id]['sequence'].encode()).hexdigest())
        sequence_info.append(nt_seqs[nt_seq_id]['sequence'])
        results[acc] = sequence_info
    return results


def main():
    """
    args[0] = str repr of path to FASTA file of query protein sequences
        (may be translated seqs from predicted ORF)
    args[1] = str repr of path to originally submitted FASTA file
        (may contain the original nucleic sequences)
    args[2] = str repr of bool, if nucleic seqs provided ('true') or not ('false')
    args[3] = str of applications
    """
    args = sys.argv[1:]
    is_fasta_check(args[0])
    sequences = get_sequences(args[0])

    if args[2] == "true":
        is_fasta_check(args[1])
        nt_seqs = get_nucleic_seqs(args[1])
        sequence_parsed = parse_nucleic(sequences, nt_seqs)
    else:
        sequence_parsed = parse(sequences, args[3])

    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
