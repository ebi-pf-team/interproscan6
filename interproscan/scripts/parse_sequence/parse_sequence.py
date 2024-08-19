import hashlib
import json
import re
import sys


NT_SEQ_ID_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=(\d+\.+\d+)\s+.+frame=\d+\s+desc=.*$")


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

def parse(sequences: dict) -> dict[str, list]:
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
    args[2] = str repr of bool, if nucleic seqs provided ('true') or not ('false')"""
    args = sys.argv[1:]
    is_fasta_check(args[0])
    sequences = get_sequences(args[0])

    if args[2] == "true":
        is_fasta_check(args[1])
        nt_seqs = get_nucleic_seqs(args[1])
        sequence_parsed = parse_nucleic(sequences, nt_seqs)
    else:
        sequence_parsed = parse(sequences)

    print(json.dumps(sequence_parsed))


if __name__ == "__main__":
    main()
