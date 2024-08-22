import sys


class NuleicError(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


def is_nucleic(seq):
    allowed_chars = ['A', 'T', 'C', 'G', 'U', 'a', 't', 'c', 'g', 'u', '*', '-']
    seq_chars = seq.upper()
    return all(char in allowed_chars for char in seq_chars)


def main(fasta):
    contains_non_nucleic_seq = 0

    with open(fasta, "r") as fh:
        current_sequence = ''
        is_seq = False

        for line in fh:
            if line.startswith('>'):
                if is_seq and not is_nucleic(current_sequence):
                    contains_non_nucleic_seq = 1
                is_seq = True
                current_sequence = ''
            else:
                current_sequence += line.strip()

        if is_seq and not is_nucleic(current_sequence):
            contains_non_nucleic_seq = 1

    if contains_non_nucleic_seq != 0:
        sys.tracebacklimit = 0
        raise NuleicError(
            """ERROR:
            The '--nucleic' flag was used, but the input FASTA file
            appears to contain at least one sequence that contains a
            non-nucleic acid residue ('A','G','C','T','*','-', case insensitive).
            Please check your input is correct.
            """
        )


if __name__ == "__main__":
    main(sys.argv[1])
