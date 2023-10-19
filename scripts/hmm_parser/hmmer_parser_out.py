import re
import argparse
import json
import members_parser.members_regex as members_regex
from cigar_alignment import cigar_alignment_parser, encode

END_OF_OUTPUT_FILE = "[ok]"
END_OF_RECORD = "//"
DOMAIN_SECTION_START = ">> "
START_OF_DOMAIN_ALIGNMENT_SECTION = "Alignments for each domain"
DOMAIN_ALIGNMENT_SECTION_START = "=="

DOMAIN_SECTION_START_PATTERN = re.compile(r"^>>\s+(\S+).*$")
DOMAIN_ALIGNMENT_LINE_PATTERN = re.compile("^\s+==\s+domain\s+(\d+)\s+.*$")
ALIGNMENT_SEQUENCE_PATTERN = re.compile("^\s+(\w+)\s+(\S+)\s+([-a-zA-Z]+)\s+(\S+)\s*$")
DOMAIN_LINE_PATTERN = re.compile(
                                "^\\s+(\\d+)\\s+[!?]\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\S+).*$")
# MODEL_ACCESSION_LINE_PATTERN = re.compile("^[^:]*:\\s+(\\w+)\\s+\\[M=(\\d+)\\].*$")


def parse(out_file):
    search_record = []
    current_sequence = None
    domain_number = None
    domains = []
    hmmer_parser_support = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    appl = out_file.split("_")[1].split(".")[0]
    member_accession = members_regex.get_accession_regex(appl)

    with open(out_file, "r") as f:
        for line in f.readlines():
            if line.startswith(END_OF_OUTPUT_FILE):
                break
            else:
                if stage == 'LOOKING_FOR_DOMAIN_DATA_LINE' and line.startswith(DOMAIN_SECTION_START):
                    stage = 'LOOKING_FOR_DOMAIN_SECTION'
                if line.startswith(END_OF_RECORD):
                    if domain_match:
                        cigar_alignment = cigar_alignment_parser(domain_match["alignment"])
                        domain_match["alignment_encoded"] = encode(cigar_alignment)
                        domains.append(domain_match)
                    sequence_match["accession"] = model_id
                    sequence_match["domains"] = domains
                    search_record.append(sequence_match)
                    hmmer_parser_support[current_sequence] = search_record
                    search_record = []
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
                        else:
                            sequence_match = get_sequence_match(line)
                            if sequence_match:
                                current_sequence = sequence_match["sequence_identifier"]
                    elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                        if line.startswith(DOMAIN_SECTION_START):
                            domain_section_header_matcher = DOMAIN_SECTION_START_PATTERN.match(line)
                            if domain_section_header_matcher:
                                domains = []
                                current_sequence = domain_section_header_matcher.group(1)
                            stage = 'LOOKING_FOR_DOMAIN_DATA_LINE'
                        # if is_tsv_pro:
                        if line.strip().startswith(DOMAIN_ALIGNMENT_SECTION_START):
                            domain_alignment_matcher = DOMAIN_ALIGNMENT_LINE_PATTERN.match(line)
                            if domain_alignment_matcher:
                                align_seq = []
                                domain_number = domain_alignment_matcher.group(1)
                        if domain_number and current_sequence:
                            alignment_sequence_pattern = ALIGNMENT_SEQUENCE_PATTERN.match(line)
                            if alignment_sequence_pattern:
                                align_seq.append(alignment_sequence_pattern.group(3))
                                domain_match["alignment"] = "".join(align_seq)

                    elif stage == 'LOOKING_FOR_DOMAIN_DATA_LINE':
                        if START_OF_DOMAIN_ALIGNMENT_SECTION in line:
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                        else:
                            match = DOMAIN_LINE_PATTERN.match(line)
                            if match:
                                domain_match = get_domain_match(match)
    return hmmer_parser_support


def get_domain_match(match):
    domain_match = {}
    domain_match["score"] = match.group(2)
    domain_match["bias"] = match.group(3)
    domain_match["cEvalue"] = match.group(4)
    domain_match["iEvalue"] = match.group(5)
    domain_match["hmm_from"] = match.group(6)
    domain_match["hmm_to"] = match.group(7)
    domain_match["hmm_bounds"] = match.group(8)
    domain_match["ali_from"] = match.group(9)
    domain_match["ali_to"] = match.group(10)
    domain_match["env_from"] = match.group(11)
    domain_match["env_to"] = match.group(12)
    domain_match["acc"] = match.group(13)
    return domain_match


def get_sequence_match(sequence_line):
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
    parser = argparse.ArgumentParser(
        description="hmmer out parser"
    )
    parser.add_argument(
        "-hmmer_file", "--hmmer_file", type=str, help="out file result of hmmer preproc")
    args = parser.parse_args()

    parse_result = parse(args.hmmer_file)
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
