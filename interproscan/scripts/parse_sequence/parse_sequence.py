import hashlib
import json
import logging
import re
import sys


NT_SEQ_ID_PATTERN = re.compile(r"^orf\d+\s+source=(.*)\s+coords=(\d+\.+\d+)\s+.+frame=\d+\s+desc=.*$")
ILLEGAL_CHARS = {
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
    "pirsr": "-",
    "prints": "-._",
    "prosite_patterns": "",
    "prosite_profiles": "-._",
    "sfld": "-._",
    "smart": "",
    "superfamily": "-",
    "signalp": "",
    "phobius": "-*._oxuzj",
}


logger = logging.getLogger()


class IllegalCharError(Exception):
    def __init__(self, fasta: str, errors: dict[str, dict[str, set]]):
        self.fasta = fasta
        self.errors = errors
        self.message = None

    def _get_error_msg(self):
        msg = f"Illegal characters detected in {self.fasta}\n"
        for seq_key in self.errors:
            for app in self.errors[seq_key]:
                app_placeholder = f" for {app}" if app != "GENERAL" else ""
                msg += (
                    f"Sequence {seq_key.split(maxsplit=1)[0]} contains illegal "
                    f"character(s){app_placeholder}: {', '.join(self.errors[seq_key][app])}\n"
                )
        self.message = msg

    def __str__(self) -> str:
        self._get_error_msg()
        return self.message


class Sequence:
    seq_key = None
    name = None
    sequence = ""

    def get_seq_key(self, line: str, nucleic: bool):
        self.seq_key = line[1:].split(maxsplit=1)[0] if nucleic else line[1:]
        self.name = line[1:] if nucleic else None

    def get_seq(self, line: str):
        self.sequence += line

    def get_error_msg(self, errors: dict[str, dict[str, str]], line: str, applications: list):
        if self.seq_key not in errors:
            errors[self.seq_key] = {}
        invalid_chars = re.findall(r"[^A-Za-z_\-\*\.]*", line)
        invalid_chars = [f"'{_char}'" for _match in set(invalid_chars)
                            for _char in _match if len(_match) > 0]
        if invalid_chars:
            if 'GENERAL' not in errors[self.seq_key]:
                errors[self.seq_key]['GENERAL'] = set()
            errors[self.seq_key]['GENERAL'] = errors[self.seq_key]['GENERAL'].union(invalid_chars)

        for app in applications:
            app_chars = [f"'{_}'" for _ in set(line).intersection(set(ILLEGAL_CHARS[app]))]
            if app_chars:
                if app.upper() not in errors[self.seq_key]:
                    errors[self.seq_key][app.upper()] = set()
                errors[self.seq_key][app.upper()] = errors[self.seq_key][app.upper()].union(app_chars)
        return errors


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
    applications: str,
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
    errors = {}
    applications = applications.split(",")
    all_illegal_chars = {chara for app in applications for chara in ILLEGAL_CHARS[app]}.union({' '})

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
                        seq_obj, sequences,
                        nucleic_seqs=nucleic_seqs, passing_nucleic=passing_nucleic
                    )

                # start with new seq
                seq_obj = Sequence()
                seq_obj.get_seq_key(line, passing_nucleic)

            else:
                if set(line.lower()).intersection(all_illegal_chars) or re.findall(r"[^A-Za-z_\-\*\.]*", line):
                    errors = seq_obj.get_error_msg(errors, line, applications)

                seq_obj.get_seq(line)

    # store the final sequence
    if seq_obj:
        sequences = store_seq(
            seq_obj, sequences,
            nucleic_seqs=nucleic_seqs, passing_nucleic=passing_nucleic
        )

    errors = {k: v for k, v in errors.items() if v != {}}
    if errors:
        sys.tracebacklimit = 0
        raise IllegalCharError(fasta_file, errors)

    return sequences
    


def main():
    """
    args[0] = str repr of path to FASTA file of query protein sequences
        (may be translated seqs from predicted ORF)
    args[1] = str repr of path to originally submitted FASTA file
        (may contain the original nucleic sequences)
    args[2] = str repr of bool, if nucleic seqs provided ('true') or not ('false')
    args[3] = str of applications
<<<<<<< HEAD
    args[4] = str repr of path for the output file
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

    with open(args[4], "w") as fh:
        json.dump(sequence_parsed, fh)
=======
    args[4] = str repr of path to write output
    """
    args = sys.argv[1:]
    nt_seqs = parse(args[1], args[3], passing_nucleic=True) if args[2] == "true" else None
    parsed_seqs = parse(args[0], args[3], nucleic_seqs=nt_seqs)
    with open(args[4], "w") as fh:
        json.dump(parsed_seqs, fh)
>>>>>>> main


if __name__ == "__main__":
    main()
