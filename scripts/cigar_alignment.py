import re


MATCH_CHAR = 'M'
INSERT_CHAR = 'I'
DELETE_CHAR = 'D'
DELETE_SYMBOL = '-'


def cigar_alignment_parser(alignment: str) -> str:
    cigar_alignment = ""
    for char in alignment:
        if char.isupper():
            cigar_alignment += MATCH_CHAR
        elif char.islower():
            cigar_alignment += INSERT_CHAR
        elif char == DELETE_SYMBOL:
            cigar_alignment += DELETE_CHAR
        else:
            raise ValueError(f"Alignment contains unrecognised character {char}")

    return cigar_alignment


def encode(cigar_alignment: str) -> str:
    encoded_alignment = ""
    prev_char = ""
    count = 0

    for char in cigar_alignment:
        if char == prev_char:
            count += 1
        else:
            if prev_char:
                encoded_alignment += str(count) + prev_char
            count = 1
            prev_char = char
    encoded_alignment += str(count) + prev_char

    return encoded_alignment


def decode(sequence: str, encoded_alignment: str) -> str:
    decoded_alignment = ""
    matches = re.findall(r'(\d+)([MID])', encoded_alignment)

    for count, operation in matches:
        count = int(count)
        if operation == MATCH_CHAR:
            decoded_alignment += sequence[:count]
            sequence = sequence[count:]
        elif operation == INSERT_CHAR:
            decoded_alignment += sequence[:count].lower()
        elif operation == DELETE_CHAR:
            decoded_alignment += DELETE_SYMBOL * count
        else:
            raise ValueError(f"Invalid CIGAR operation: {operation}")

    return decoded_alignment
