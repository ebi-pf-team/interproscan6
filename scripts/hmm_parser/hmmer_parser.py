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


def parse_out(out_file):
    # TODO: this parse is based on i5 parser, maybe it can be simplified (similar to domtbl parser)
    search_record = {}
    current_domain = None
    domains = {}
    hmmer3ParserSupport = []
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    appl = out_file.split("_")[1].split(".")[0]
    member_accession = members_regex.get_accession_regex(appl)

    with open(out_file, "r") as f:
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
                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.strip() == "":
                            if search_record:
                                stage = 'LOOKING_FOR_DOMAIN_SECTION'
                            else:
                                stage = 'FINISHED_SEARCHING_RECORD'
                            current_domain = None
                        else:
                            sequence_match = get_sequence_match(line)
                            if sequence_match:
                                current_sequence = sequence_match["sequenceIdentifier"]
                                try:
                                    search_record[current_sequence].append(sequence_match)
                                except:
                                    search_record[current_sequence] = [sequence_match]
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


def parse_domtbl(domtbl_file):
    sequence_matches = {}
    with open(domtbl_file, "r") as f:
        current_seq = None
        acc = []
        domains = []
        for line in f.readlines():
            if not line.startswith(COMMENT_LINE):
                info = line.split()
                if info[0] != current_seq:
                    if current_seq:
                        sequence_matches[current_seq] = {
                            "accession": info[1],
                            "tlen": int(info[2]),
                            "acc_matches": acc
                        }
                    acc = []
                    current_seq = info[0]
                if info[9] == info[10]:
                    domains.append(get_domain(info))
                    acc.append(get_accession(info, domains))
                    domains = []
                else:
                    domains.append(get_domain(info))

            if current_seq:
                sequence_matches[current_seq] = {
                    "accession": info[1],
                    "tlen": int(info[2]),
                    "acc_matches": acc
                }
    return sequence_matches


def get_accession(info, domains):
    acc_info = {
        "query_name": info[3],
        "accession": info[4],
        "qlen": int(info[5]),
        "e_value": float(info[6]),
        "score": float(info[7]),
        "bias": float(info[8]),
        "domains": domains
    }
    return acc_info


def get_domain(info):
    domain_info = {
        "cEvalue": info[11],
        "iEvalue": info[12],
        "score": info[13],
        "bias": info[14],
        "hmm_from": info[15],
        "hmm_to": info[16],
        "ali_from": info[17],
        "ali_to": info[18],
        "env_from": info[19],
        "env_to": info[20],
        "acc": info[21],
        "description_of_target": info[22]
    }
    return domain_info


def main():
    parser = argparse.ArgumentParser(
        description="hmmer parser"
    )

    parser.add_argument(
        "-out", "--preproc_out", type=str, help="out file result of hmmer preproc")
    parser.add_argument(
        "-domtbl", "--preproc_domtbl", type=str, help="dtbl file result of hmmer preproc")
    args = parser.parse_args()

    parse_dtbl_result = parse_domtbl(args.preproc_domtbl)
    print(json.dumps(parse_dtbl_result, indent=2))


if __name__ == "__main__":
    main()
