import hashlib
import json
import re
import sys


NT_SEQ_ID_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=(\d+\.+\d+)\s+.+frame=\d+\s+desc=.*$")


class Sequence:
    seq_key = None
    name = None
    sequence = ""

    def get_seq_key(self, line: str, nucleic: bool):
        self.seq_key = line[1:].split(maxsplit=1)[0] if nucleic else line[1:]
        self.name = line[1:] if nucleic else None

    def get_seq(self, line: str):
        self.sequence += line


def store_seq(
    seq_obj: Sequence,
    sequences: dict[str, dict[str, Sequence]],
    nucleic_seqs=None,
    passing_nucleic=False
):
    acc = seq_obj.seq_key.split(maxsplit=1)[0]
    if passing_nucleic:
        sequences[acc] = seq_obj
    else:
        sequences[acc] = {
            'seq_id': seq_obj.seq_key,
            'name': seq_obj.name,
            'sequence': seq_obj.sequence,
            'md5': hashlib.md5(seq_obj.sequence.encode()).hexdigest(),
            'length': len(seq_obj.sequence)
        }
        if nucleic_seqs:
            nt_acc = NT_SEQ_ID_PATTERN.match(seq_obj.seq_key)
            nt_acc = nt_acc.group(1)
            # Passing ORF FASTA after input nucleic FASTA
            # Add additional nt data for mapping ORF to the nt seq in the output
            sequences[acc]['nt_seq_id'] = nucleic_seqs[nt_acc].seq_key
            sequences[acc]['nt_name'] = nucleic_seqs[nt_acc].name
            sequences[acc]['nt_sequence'] = nucleic_seqs[nt_acc].sequence
            sequences[acc]['nt_md5'] = hashlib.md5(
                nucleic_seqs[nt_acc].sequence.encode()).hexdigest()
    return sequences


def parse(
    fasta_file: str,
    passing_nucleic=False,
    nucleic_seqs=None
) -> dict[str, dict[str, any]]:
    """
    :param fasta_file: str repr of path to FASTA file
    :param passing_nucleic: bool, passing nucleic acid seq FASTA
    :param applications: str repr of list of applications to be used in IPS6 run
    :param nucleic_seqs: dict of seqs from input nucleic acid seq FASTA
    """
    seq_obj = None
    sequences = {}

    with open(fasta_file, "r") as fh:
        if not fh.readline().strip().startswith(">"):
            raise ValueError(f"{fasta_file} is not in FASTA format")

    with open(fasta_file, "r") as fh:
        for line in fh:
            line = line.strip()

            if line.startswith(">"):
                # store completed seq
                if seq_obj:
                    sequences = store_seq(
                        seq_obj,
                        sequences,
                        nucleic_seqs=nucleic_seqs,
                        passing_nucleic=passing_nucleic
                    )

                # start with new seq
                seq_obj = Sequence()
                seq_obj.get_seq_key(line, passing_nucleic)

            else:
                seq_obj.get_seq(line)

    # store the final sequence
    if seq_obj:
        sequences = store_seq(
            seq_obj,
            sequences,
            nucleic_seqs=nucleic_seqs,
            passing_nucleic=passing_nucleic
        )

    return sequences


def main():
    """
    args[0] = str repr of path to FASTA file of query protein sequences
        (may be translated seqs from predicted ORF)
    args[1] = str repr of path to originally submitted FASTA file
        (may contain the original nucleic sequences)
    args[2] = str repr of bool, if nucleic seqs provided ('true') or not ('false')
    args[3] = str repr of path to write output
    """
    args = sys.argv[1:]
    nt_seqs = parse(args[1], passing_nucleic=True) if args[2] == "true" else None
    parsed_seqs = parse(args[0], nucleic_seqs=nt_seqs)
    with open(args[3], "w") as fh:
        json.dump(parsed_seqs, fh)


if __name__ == "__main__":
    main()
