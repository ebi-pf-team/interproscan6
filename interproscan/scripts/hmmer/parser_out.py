import json
import re
import sys

DOMAIN_SECTION_START_PATTERN = re.compile(r"^>>\s+(\S+).*$")
DOMAIN_ALIGNMENT_LINE_PATTERN = re.compile(r"^\s+==\s+domain\s+(\d+)\s+.*$")
ALIGNMENT_SEQUENCE_PATTERN = re.compile(r"^\s+(\S+)\s+(\S+)\s+([-a-zA-Z]+)\s+(\S+)\s*$")  # replacing (\w+) with (\S+) and adding if to ignore current sequence
HMMALIGN_SEQUENCE_PATTERN = re.compile(r"^\s+(\S+)\s+(\S+)\s+([-a-zA-Z.]+)\s+(\S+)\s*$")
DOMAIN_LINE_PATTERN = re.compile(
    "^\\s+(\\d+)\\s+[!?]\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\S+).*$"
)

MATCH_CHAR = 'M'
INSERT_CHAR = 'I'
DELETE_CHAR = 'D'
DELETE_SYMBOL = '-'


def get_accession_regex(appl: str) -> re.Pattern:
    if appl.upper() == "ANTIFAM":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(ANF\d{5})\s*$")
    if appl.upper() == "FUNFAM":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+([\w\d\.\-]+)\s*\[M=\d+\]$")
    if appl.upper() == "GENE3D":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+([\w\d\.\-_]+)\s*\[M=\d+\]$")
    if appl.upper() == "HAMAP":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(MF_\d+)\s+(\[M=\d+\])$")
    if appl.upper() == "NCBIFAM":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+((TIGR|NF)\d+)\.\d+$")
    if appl.upper() == "PANTHER":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(PTHR[^\s]+)\s*\[M=\d+\]$")
    if appl.upper() == "PFAM":
        return re.compile(r"(^[^:]*:)\s+(\w+).*$")
    if appl.upper() == "PIRSF":
        return re.compile(r"^(Accession:)\s+(PIRSF\d+)$")
    if appl.upper() == "PIRSR":
        return re.compile(r"^(Accession:)\s+(PIRSR\d+-\d+)$")
    if appl.upper() == "SFLD":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(SFLD[^\s]+)\s*$")


def parse(out_file: str) -> dict:
    """Main parser for the hmmer.out file.

    :parama out_file: str representation of path to HMMER.out file
    """
    member_db = out_file.split("/")[-1].split("._.")[0]
    current_sequence = None
    current_domain = None
    sequence_match = {}
    domain_match = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    member_accession = get_accession_regex(member_db)
    model_id = ""

    with open(out_file, "r") as f:
        for line in f.readlines():
            if line.startswith("[ok]"):
                break
            if line.startswith("#"):
                pass
            else:
                if stage == 'LOOKING_FOR_DOMAIN_DATA_LINE' and line.startswith(">> "):
                    stage = 'LOOKING_FOR_DOMAIN_SECTION'
                if line.startswith("//"):
                    if domain_match:
                        sequence_match[current_sequence][model_id] = add_domain_match(
                            domain_match,
                            sequence_match[current_sequence][model_id]
                        )
                        domain_match = {}
                    stage = "LOOKING_FOR_METHOD_ACCESSION"
                else:
                    if stage == 'LOOKING_FOR_METHOD_ACCESSION':
                        if line.startswith("Accession:") or \
                            line.startswith("Query sequence:") or \
                            (line.startswith("Query:") and member_db.lower() in ["funfam", "gene3d", "hamap", "panther"]):
                            stage = 'LOOKING_FOR_SEQUENCE_MATCHES'
                            model_ident_pattern = member_accession.match(line)

                            if model_ident_pattern:
                                model_id = model_ident_pattern.group(2) if \
                                    member_db != "panther" else \
                                    model_ident_pattern.group(2).replace(".orig.30.pir", "")

                        if line.startswith("Query:"):
                            try:
                                qlen = line.split("[")[1].split("]")[0].replace("M=", "")
                            except IndexError:  # e.g. line = "Query:       5xqwL01-i2]"
                                qlen = ""

                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.strip() == "":
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                            current_domain = None
                            current_sequence = None
                        else:
                            sequence_match = get_sequence_match(
                                sequence_match,
                                line,
                                model_id,
                                member_db,
                                qlen
                            )

                    elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                        if line.startswith(">> "):
                            domain_section_header_matcher = DOMAIN_SECTION_START_PATTERN.match(line)
                            if domain_section_header_matcher:
                                if domain_match:
                                    sequence_match[current_sequence][model_id] = add_domain_match(
                                        domain_match,
                                        sequence_match[current_sequence][model_id]
                                    )
                                    domain_match = {}
                                current_sequence = domain_section_header_matcher.group(1)
                            stage = 'LOOKING_FOR_DOMAIN_DATA_LINE'
                        if line.strip().startswith("=="):
                            domain_alignment_matcher = DOMAIN_ALIGNMENT_LINE_PATTERN.match(line)
                            if domain_alignment_matcher:
                                align_seq = []
                                hmm_seq = []
                                current_domain = domain_alignment_matcher.group(1)
                        if current_domain and current_sequence:
                            alignment_sequence_pattern = ALIGNMENT_SEQUENCE_PATTERN.match(line)
                            if alignment_sequence_pattern:
                                if alignment_sequence_pattern.group(1) == current_sequence:
                                    align_seq.append(alignment_sequence_pattern.group(3))
                                    domain_match[current_domain]["alignment"] = "".join(align_seq)
                            if member_db.upper() == "PIRSR":
                                hmmalign_sequence_pattern = HMMALIGN_SEQUENCE_PATTERN.match(line)
                                if hmmalign_sequence_pattern:
                                    if hmmalign_sequence_pattern.group(1) != current_sequence:
                                        hmm_seq.append(hmmalign_sequence_pattern.group(3))
                                        domain_match[current_domain]["hmm_alignment"] = "".join(hmm_seq)
                    elif stage == 'LOOKING_FOR_DOMAIN_DATA_LINE':
                        if "Alignments for each domain" in line:
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                        else:
                            match = DOMAIN_LINE_PATTERN.match(line)
                            if match:
                                domain_number = match.group(1)
                                try:
                                    domain_match[domain_number] = get_domain_match(match, member_db, qlen)
                                except KeyError:
                                    domain_match = {domain_number: get_domain_match(match, member_db, qlen)}

    return sequence_match


def add_domain_match(domain_match: dict, current_match: dict) -> dict:
    """
    Add domain match as a match 'location'.

    :param domain_match: dict representing a domain hit from HMMER
    :param current_match: dict representing current signature match
    """
    for domain_key, domain_value in domain_match.items():
        cigar_alignment = cigar_alignment_parser(domain_match[domain_key]["alignment"])
        domain_match[domain_key]["cigar_alignment"] = encode(cigar_alignment)
        current_match["locations"].append(domain_match[domain_key])
    return current_match


def get_domain_match(match: re.Match, member_db: str, qlen: str) -> dict:
    """
    Retrieve data about a domain match

    :param match: re group object matching a line containing a domain hit
    :param member_db: name of the member database
    :param qlen: query length
    """
    domain_match = {}

    hmm_bound_pattern = {
        "[]": "COMPLETE",
        "[.": "N_TERMINAL_COMPLETE",
        ".]": "C_TERMINAL_COMPLETE",
        "..": "INCOMPLETE"
    }
    domain_match["start"] = int(match.group(9)) if member_db.upper() != "GENE3D" else int(match.group(11))  # ali coord from
    domain_match["end"] = int(match.group(10)) if member_db.upper() != "GENE3D" else int(match.group(12))  # ali coord to
    domain_match["representative"] = ""
    domain_match["hmmStart"] = int(match.group(6))  # hmm coord from
    domain_match["hmmEnd"] = int(match.group(7))   # hmm coord to
    domain_match["hmmLength"] = int(qlen)  # qlen
    domain_match["rawHmmBounds"] = match.group(8)
    domain_match["hmmBounds"] = hmm_bound_pattern[match.group(8)]
    domain_match["evalue"] = float(match.group(5))  # Independent e-value
    domain_match["score"] = float(match.group(2))   # bit score
    domain_match["envelopeStart"] = match.group(11)    # env coord from
    domain_match["envelopeEnd"] = match.group(12)  # env coord to

    return domain_match


def get_sequence_match(
    sequence_match: dict,
    sequence_line: str,
    model_id: str,
    member_db: str,
    qlen: str
) -> dict:
    """Retrieve data about the query sequence"""
    SEQUENCE_LINE_PATTERN = re.compile(
        "^\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+(\\S+).*$"
    )
    match = SEQUENCE_LINE_PATTERN.match(sequence_line)
    if match:
        match_info = {
            "accession": model_id,
            "evalue": float(match.group(1)),
            "score": float(match.group(2)),
            "qlen": int(qlen),
            "bias": match.group(3),
            "member_db": member_db,
            "model-ac": model_id.split(":")[0].split(".")[0],
            "locations": []
        }
        sequence_id = match.group(4)
        if sequence_id not in sequence_match:
            sequence_match[sequence_id] = {}

        if model_id not in sequence_match[sequence_id]:
            sequence_match[sequence_id][model_id] = {}

        sequence_match[sequence_id][model_id].update(match_info)

    return sequence_match


def cigar_alignment_parser(alignment: str) -> str:
    """Convert alignment from HMMER3 and HMMER2 to cigar alignment

    :param alignment: The query protein sequence from the alignment
    """
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


def main():
    """
    :args 0: str repr of path to hmmer file to be parsed
    :args 1: str repr of path to write output
    """
    args = sys.argv[1:]
    parse_result = parse(args[0])
    output = args[1]

    with open(output, "w") as fh:
        json.dump(parse_result, fh)


if __name__ == "__main__":
    main()
