import re
import argparse
import json
import members_parser.members_regex as members_regex

COMMENT_LINE = "#"
END_OF_OUTPUT_FILE = "[ok]"
END_OF_RECORD = "//"
DOMAIN_SECTION_START = ">> "
START_OF_DOMAIN_ALIGNMENT_SECTION = "Alignments for each domain"
DOMAIN_ALIGNMENT_SECTION_START = "=="

DOMAIN_SECTION_START_PATTERN = re.compile(r"^>>\s+(\S+).*$")
DOMAIN_ALIGNMENT_LINE_PATTERN = re.compile("^\\s+==\\s+domain\\s+(\\d+)\\s+.*$")
ALIGNMENT_SEQUENCE_PATTERN = re.compile("^\\s+(\\w+)\\s+(\\S+)\\s+([-a-zA-Z]+)\\s+(\\S+)\\s*$")
MODEL_ACCESSION_LINE_PATTERN = re.compile("^[^:]*:\\s+(\\w+)\\s+\\[M=(\\d+)\\].*$")


def parse(rpsbproc_file, appl):
    search_record = {}
    current_domain = None
    domains = {}
    sequences = {}
    hmmer3ParserSupport = []
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    member_accession = members_regex.get_accession_regex(appl)

    with open(rpsbproc_file, "r") as f:
        for line in f.readlines():
            if line.startswith(COMMENT_LINE):
                pass
            elif line.startswith(END_OF_OUTPUT_FILE):
                break
            else:
                if stage == 'LOOKING_FOR_DOMAIN_DATA_LINE' and line.startswith(DOMAIN_SECTION_START):
                    stage = 'LOOKING_FOR_DOMAIN_SECTION'
                if line.startswith(END_OF_RECORD):
                    if search_record:
                        hmmer3ParserSupport.append(search_record)
                    search_record = {}
                    stage = "LOOKING_FOR_METHOD_ACCESSION"
                else:
                    if stage == 'LOOKING_FOR_METHOD_ACCESSION':
                        if line.startswith("Accession:") or line.startswith("Query sequence:"):
                            stage = 'LOOKING_FOR_SEQUENCE_MATCHES'
                            model_ident_pattern = member_accession.match(line)
                            if model_ident_pattern:
                                model_id = model_ident_pattern.group(1)
                                search_record["signature_acc"] = model_id
                                # model_length = model_ident_pattern.group(2)
                                # search_record["model_length"] = model_length
                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.strip() == "":
                            if search_record:
                                stage = 'LOOKING_FOR_DOMAIN_SECTION'
                            else:
                                stage = 'FINISHED_SEARCHING_RECORD'
                            current_domain = None
                            current_sequence = None
                        else:
                            sequence_match = get_sequence_match(line)
                            if sequence_match:
                                current_sequence = sequence_match["sequenceIdentifier"]
                                try:
                                    sequences[current_sequence].append(sequence_match)
                                except:
                                    sequences[current_sequence] = [sequence_match]
                            search_record["sequence_match"] = sequences
                    elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                        if line.startswith(DOMAIN_SECTION_START):
                            domains = {}
                            match = DOMAIN_SECTION_START_PATTERN.match(line)
                            if match:
                                current_domain = match.group(1)
                            stage = 'LOOKING_FOR_DOMAIN_DATA_LINE'
                    elif stage == 'LOOKING_FOR_DOMAIN_DATA_LINE':
                        if line == START_OF_DOMAIN_ALIGNMENT_SECTION:
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                        else:
                            domain_match = get_domain_match(line)
                            current_domain = match.group(1)

                        if domain_match:
                            try:
                                domains[current_domain].append(domain_match)
                            except:
                                domains[current_domain] = [domain_match]
                        search_record["domain_match"] = domains
    return hmmer3ParserSupport


def get_domain_match(domain_line):
    domain_match = {}
    DOMAIN_LINE_PATTERN = re.compile(
        "^\\s+(\\d+)\\s+[!?]\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\d+)\\s+(\\d+)\\s+\\S+\\s+(\\S+).*$")
    match = DOMAIN_LINE_PATTERN.match(domain_line)
    if match:
        domain_match["score"] = match.group(2)
        domain_match["bias"] = match.group(3)
        domain_match["cEvalue"] = match.group(4)
        domain_match["iEvalue"] = match.group(5)
        domain_match["hmmfrom"] = match.group(6)
        domain_match["hmmto"] = match.group(7)
        domain_match["hmmBounds"] = match.group(8)
        domain_match["aliFrom"] = match.group(9)
        domain_match["aliTo"] = match.group(10)
        domain_match["envFrom"] = match.group(11)
        domain_match["envTo"] = match.group(12)
        domain_match["acc"] = match.group(13)
    return domain_match


def get_sequence_match(sequence_line):
    sequence_match = {}
    SEQUENCE_LINE_PATTERN = re.compile("^\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+(\\S+).*$")
    match = SEQUENCE_LINE_PATTERN.match(sequence_line)
    if match:
        sequence_match["eValue"] = match.group(1)
        sequence_match["score"] = match.group(2)
        sequence_match["bias"] = match.group(3)
        sequence_match["sequenceIdentifier"] = match.group(4)
    return sequence_match


def main():
    parser = argparse.ArgumentParser(
        description="hmmer parser"
    )
    parser.add_argument(
        "-seq", "--sequences", type=str, help="fasta file with sequences"
    )
    parser.add_argument(
        "-preproc", "--preproc", type=str, help="file result of hmmer preproc")
    parser.add_argument(
        "-appl", "--application", type=str, help="name of member database")
    args = parser.parse_args()

    hmmer_parse_result = parse(args.preproc, args.application)
    print(json.dumps(hmmer_parse_result))


if __name__ == "__main__":
    main()
