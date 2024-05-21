import json
import re
import sys
from cigar_alignment import cigar_alignment_parser, encode


DOMAIN_SECTION_START_PATTERN = re.compile(r"^>>\s+(\S+).*$")
DOMAIN_ALIGNMENT_LINE_PATTERN = re.compile(r"^\s+==\s+domain\s+(\d+)\s+.*$")
ALIGNMENT_SEQUENCE_PATTERN = re.compile(r"^\s+(\w+)\s+(\S+)\s+([-a-zA-Z]+)\s+(\S+)\s*$")
DOMAIN_LINE_PATTERN = re.compile(
                                "^\\s+(\\d+)\\s+[!?]\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\S+).*$")


def get_accession_regex(appl: str) -> re.Pattern:
    if appl.upper() == "ANTIFAM":
        return re.compile(r"^Accession:\s+(ANF\d{5})\s*$")
    if appl.upper() == "NCBIFAM":
        return re.compile(r"^Accession:\s+((TIGR|NF)\d+)\.\d+$")
    if appl.upper() == "SFLD":
        return re.compile(r"^Accession:\s+(SFLD[^\s]+)\s*$")


def parse(out_file: str) -> dict:
    current_sequence = None
    current_domain = None
    domain_match = {}
    hmmer_parser_support = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    appl = out_file.split("_")[1].split(".")[0]
    member_accession = get_accession_regex(appl)

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
                        for domain_key, domain_value in domain_match.items():
                            cigar_alignment = cigar_alignment_parser(domain_match[domain_key]["alignment"])
                            domain_match[domain_key]["cigar_alignment"] = encode(cigar_alignment)
                        sequence_match["sequence"] = current_sequence
                        sequence_match["domains"] = domain_match
                        domain_match = {}
                        if current_sequence in hmmer_parser_support:
                            hmmer_parser_support[current_sequence].append(sequence_match)
                        else:
                            hmmer_parser_support[current_sequence] = [sequence_match]
                    sequence_match = {}

                    stage = "LOOKING_FOR_METHOD_ACCESSION"
                else:
                    if stage == 'LOOKING_FOR_METHOD_ACCESSION':
                        if line.startswith("Accession:") or line.startswith("Query sequence:"):
                            stage = 'LOOKING_FOR_SEQUENCE_MATCHES'
                            model_ident_pattern = member_accession.match(line)
                            if model_ident_pattern:
                                model_id = model_ident_pattern.group(1)
                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.strip() == "":
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                            current_domain = None
                            current_sequence = None
                        else:
                            sequence_match = get_sequence_match(line)
                            if sequence_match:
                                current_sequence = sequence_match["sequence_identifier"]
                    elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                        if line.startswith(">> "):
                            domain_section_header_matcher = DOMAIN_SECTION_START_PATTERN.match(line)
                            if domain_section_header_matcher:
                                current_sequence = domain_section_header_matcher.group(1)
                            stage = 'LOOKING_FOR_DOMAIN_DATA_LINE'
                        if line.strip().startswith("=="):
                            domain_alignment_matcher = DOMAIN_ALIGNMENT_LINE_PATTERN.match(line)
                            if domain_alignment_matcher:
                                align_seq = []
                                current_domain = domain_alignment_matcher.group(1)
                        if current_domain and current_sequence:
                            alignment_sequence_pattern = ALIGNMENT_SEQUENCE_PATTERN.match(line)
                            if alignment_sequence_pattern:
                                align_seq.append(alignment_sequence_pattern.group(3))
                                domain_match[current_domain]["alignment"] = "".join(align_seq)

                    elif stage == 'LOOKING_FOR_DOMAIN_DATA_LINE':
                        if "Alignments for each domain" in line:
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                        else:
                            match = DOMAIN_LINE_PATTERN.match(line)
                            if match:
                                domain_number = match.group(1)
                                domain_match[domain_number] = get_domain_match(match)
    return hmmer_parser_support


def get_domain_match(match: re.Match) -> dict:
    domain_match = {}
    hmm_bound_pattern = {"[]": "Complete", "[.": "N-terminal complete", ".]": "C-terminal complete", "..": "Incomplete"}
    domain_match["score"] = match.group(2)
    domain_match["bias"] = match.group(3)
    domain_match["cEvalue"] = match.group(4)
    domain_match["iEvalue"] = match.group(5)
    domain_match["hmm_from"] = match.group(6)
    domain_match["hmm_to"] = match.group(7)
    domain_match["hmm_bounds"] = match.group(8)
    domain_match["hmm_bounds_parsed"] = hmm_bound_pattern[match.group(8)]
    domain_match["ali_from"] = match.group(9)
    domain_match["ali_to"] = match.group(10)
    domain_match["env_from"] = match.group(11)
    domain_match["env_to"] = match.group(12)
    domain_match["accession"] = match.group(13)
    return domain_match


def get_sequence_match(sequence_line: str) -> dict:
    sequence_match = {}
    SEQUENCE_LINE_PATTERN = re.compile("^\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+(\\S+).*$")
    match = SEQUENCE_LINE_PATTERN.match(sequence_line)
    if match:
        sequence_match["e_value"] = match.group(1)
        sequence_match["score"] = match.group(2)
        sequence_match["bias"] = match.group(3)
        sequence_match["sequence_identifier"] = match.group(4)
    return sequence_match


def main():
    """
    :args 0: str repr of path to hmmer file to be parsed
    """
    args = sys.argv[1:]
    parse_result = parse(args[0])

    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
