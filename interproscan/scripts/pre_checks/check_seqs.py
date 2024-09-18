import re
import sys

from pathlib import Path


LEGAL_NUCLEIC_CHARS = ['A', 'T', 'C', 'G', 'U', 'a', 't', 'c', 'g', 'u', '*', '-']
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
    "signalp_euk": "",
    "phobius": "-*._oxuzj",
}


class raise InvalidInputError(Exception):
    def __init__(self, outpath: Path, message: str, errors: dict[str, dict|set]):
        self.outpath = outpath
        self.message = message
        self.errors = errors

    def _write_errors(self):
        """Write out file listing all character violations"""
        self.outpath.parent.mkdir(parents=True, exist_ok=True)

        error_txt = ""
        if self.errors['non-nucleic-chars']:
            error_txt += "Non-nucleic acid residues were found in the following sequences:\n"
            for seq_key in self.errors['non-nucleic-chars']:
                error_txt += f"{seq_key}\n"
            error_txt += "\n"
        if self.errors['illegal-chars']:
            error_txt += "Illegal characters were detected in the following sequences:\n"
            for seq_key, seq_errors in self.errors['illegal-chars'].items():
                for app in seq_errors:
                    app_placeholder = f" for {app}" if app != "GENERAL" else ""
                    error_txt += (
                        f"Sequence {seq_key.split(maxsplit=1)[0]} contains illegal "
                        f"character(s){app_placeholder}: {', '.join(seq_errors[app])}\n"
                    )

        with open(self.outpath, "w") as fh:
            fh.write(error_txt)

    def __str__(self):
        self._write_errors()
        return self.message


class Sequence():
    seq_key = None
    sequence = ""

    def get_seq_key(self, line: str):
        self.seq_key = line.strip(">")

    def get_seq(self, line: str):
        self.sequence += line

    def is_nucleic(self):
        "Check for non-nucleic chars"
        seq_chars = set(self.sequence.upper())
        return all(char in LEGAL_NUCLEIC_CHARS for char in seq_chars)

    def check_illegal_chars(
        self,
        applications: list,
        passing_nucleic: bool
    ):
        errors = {}

        if passing_nucleic:
            invalid_char_pattern = r"[^A-Za-z_\.]*"
        else:
            invalid_char_pattern = r"[^A-Za-z_\-\*\.]*"

        invalid_chars = re.findall(invalid_char_pattern, self.sequence)
        invalid_chars = [f"'{_char}'" for _match in set(invalid_chars)
                            for _char in _match if len(_match) > 0]
        if invalid_chars:
            if 'GENERAL' not in errors:
                errors['GENERAL'] = set()
            errors['GENERAL'] = errors['GENERAL'].union(invalid_chars)

        for app in applications:
            app_chars = [f"'{_}'" for _ in set(self.sequence).intersection(set(ILLEGAL_CHARS[app]))]
            if app_chars:
                # 'u; is an illegal char for phobius but allow it when
                # passing a nucleic acid sequence
                if passing_nucleic and app.lower() == 'phobius' and "'u'" in app_chars:
                    app_chars.remove("'u'")
                    continue
                if self.seq_key not in errors:
                    errors = {}
                if app.upper() not in errors:
                    errors[app.upper()] = set()
                errors[app.upper()] = errors[app.upper()].union(app_chars)

        return errors


def main(fasta: str, applications: str, nucleic: bool, outdir: str):
    """Test if the fasta file contains illegal chars

    :param fasta: str repr of path to input fasta file
    :param applications: str repr of list of applications
    :param nucleic: bool, should the fasta file contain nucleic seqs
    """
    applications = set(applications.split(","))
    errors = {
        'non-nucleic-chars': set(),
        'illegal-chars': {},
    }
    outpath = Path(outdir) / "illegal_characters.err"

    with open(fasta, "r") as fh:
        current_seq = None

        for line in fh:
            if line.startswith('>'):
                if current_seq:
                    # parse the previous sequence

                    # check if nucleic seq contains non-nucleic residues
                    if nucleic:
                        not_nucleic = current_seq.is_nucleic()
                        if not_nucleic:
                            errors['non-nucleic-chars'].add(current_seq.seq_key)

                    # check for general and application specific illegal chars
                    illegal_chars_detected = current_seq.check_illegal_chars(applications, nucleic)
                    if illegal_chars_detected:
                        if current_seq.seq_key not in errors:
                            errors['illegal-chars'][current_seq.seq_key] = {}
                        errors['illegal-chars'][current_seq.seq_key] = illegal_chars_detected

                current_seq = Sequence()
                current_seq.get_seq_key(line)

            else:
                current_seq.get_seq(line)

    # pass the last sequence of the file

    # check if nucleic seq contains non-nucleic residues
    if nucleic:
        not_nucleic = current_seq.is_nucleic()
        if not_nucleic:
            errors['non-nucleic-chars'].add(current_seq.seq_key)

    # check for general and application specific illegal chars
    illegal_chars_detected = current_seq.check_illegal_chars(applications, nucleic)
    if illegal_chars_detected:
        if current_seq.seq_key not in errors:
            errors['illegal-chars'][current_seq.seq_key] = {}
        errors['illegal-chars'][current_seq.seq_key] = illegal_chars_detected

    if errors['non-nucleic-chars'] and errors['illegal-chars']:
        # both detected
        message = (
            "::ERROR::\n"
            f"Illegal characters were detected in {fasta}\n"
            "and the '--nucleic' flag was used, but the input FASTA file\n"
            "appears to contain non-nucleic acid residues ('A','G','C','T','*','-', case insensitive)\n"
            f"All violations are listed in {outpath}\n"
            "Please check your input is correct."
        )
        raise InvalidInputError(
            outpath,
            message,
            errors
        )
    elif errors['non-nucleic-chars'] and not errors['illegal-chars']:
        # only non-nuleic chars detects
        message = (
            "::ERROR::\n"
            f"The '--nucleic' flag was used, but the input {fasta} file\n"
            "appears to contain non-nucleic acid residues ('A','G','C','T','*','-', case insensitive)\n"
            f"All violations are listed in {outpath}\n"
            "Please check your input is correct."
        )
        raise InvalidInputError(
            outpath,
            message,
            errors
        )
    elif not errors['non-nucleic-chars'] and errors['illegal-chars']:
        # only illegal chars detects
        message = (
            "::ERROR::\n"
            f"Illegal characters were detected in {fasta}\n"
            f"All violations are listed in {outpath}\n"
            "Please check your input is correct."
        )
        raise InvalidInputError(
            outpath,
            message,
            errors
        )


if __name__ == "__main__":
    main(sys.argv[1])
