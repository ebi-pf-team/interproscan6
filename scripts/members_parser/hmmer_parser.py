import re

COMMENT_LINE = "#"
END_OF_OUTPUT_FILE = "[ok]"

END_OF_RECORD = "//"
DOMAIN_SECTION_START = ">> "
START_OF_DOMAIN_ALIGNMENT_SECTION = "Alignments for each domain"
DOMAIN_ALIGNMENT_SECTION_START = "  =="

DOMAIN_SECTION_START_PATTERN = re.compile(r"^>>\s+(\S+).*$")
DOMAIN_ALIGNMENT_LINE_PATTERN = re.compile("^\\s+==\\s+domain\\s+(\\d+)\\s+.*$")
ALIGNMENT_SEQUENCE_PATTERN = re.compile("^\\s+(\\w+)\\s+(\\S+)\\s+([-a-zA-Z]+)\\s+(\\S+)\\s*$")


def parse(rpsbproc_file, accession_re):
    search_record = {}
    current_sequence_identifier = None
    domains = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    hmmer3ParserSupport = []

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
                        if line.startswith(DOMAIN_SECTION_START):
                            stage = 'LOOKING_FOR_SEQUENCE_MATCHES'
                            match = DOMAIN_SECTION_START_PATTERN.match(line)
                            if match:
                                current_sequence_identifier = match.group(1)
                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.strip() == "":
                            if search_record:
                                stage = 'LOOKING_FOR_DOMAIN_SECTION'
                                current_sequence_identifier = None
                            else:
                                stage = 'FINISHED_SEARCHING_RECORD'
                        else:
                            domain_match = get_domain_match(line)
                            if domain_match:
                                try:
                                    search_record[current_sequence_identifier].append(domain_match)
                                except:
                                    search_record[current_sequence_identifier] = [domain_match]
                    # elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                    #     if line.startswith(DOMAIN_SECTION_START):
                    #         match = DOMAIN_SECTION_START_PATTERN.match(line)
                    #         if match:
                    #             domains.clear()
                    #             current_sequence_identifier = match.group(1)
                    #             stage = 'LOOKING_FOR_DOMAIN_DATA_LINE'
                    # elif stage == 'LOOKING_FOR_DOMAIN_DATA_LINE':
                    #     match = get_domain_match(line)
                        # if line.strip() == START_OF_DOMAIN_ALIGNMENT_SECTION:
                        #     stage = 'LOOKING_FOR_DOMAIN_SECTION'
                        # elif match:
                        #     domain_match = get_domain_match(match)
                        #     domain_line_id = match.group(1)
                        #     search_record["domain_match"] = (current_sequence_identifier, domain_match)
                        #     domains[domain_line_id] = domain_match
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


if __name__ == "__main__":  # this main now is just to test directly, on the flow of nextflow we don't need
    result = parse("./work/b7/cf22b37d1c9adfc5c676dabb3c8142/hmmer_ncbifam.hmm_best_to_test.1.fasta.out")
    print(result)
#     the output now is like this:
#       [{'Q97R95': [{'score': '446.0', 'bias': '3.0', 'cEvalue': '5.5e-138', 'iEvalue': '8.5e-134', 'hmmfrom': '1', 'hmmto': '359', 'hmmBounds': '[.', 'aliFrom': '4', 'aliTo': '361', 'envFrom': '4', 'envTo': '365', 'acc': '0.98'}]},
#       {'UPI00043D6473': [{'score': '773.5', 'bias': '0.0', 'cEvalue': '1.7e-236', 'iEvalue': '1.3e-232', 'hmmfrom': '24', 'hmmto': '597', 'hmmBounds': '..', 'aliFrom': '100', 'aliTo': '670', 'envFrom': '83', 'envTo': '671', 'acc': '0.92'},
#                           {'score': '662.3', 'bias': '7.9', 'cEvalue': '8.1e-203', 'iEvalue': '6.3e-199', 'hmmfrom': '141', 'hmmto': '597', 'hmmBounds': '..', 'aliFrom': '723', 'aliTo': '1164', 'envFrom': '709', 'envTo': '1165', 'acc': '0.96'}]},
#       {'A0B6J9': [{'score': '370.1', 'bias': '2.4', 'cEvalue': '3.3e-115', 'iEvalue': '5.1e-111', 'hmmfrom': '1', 'hmmto': '253', 'hmmBounds': '[.', 'aliFrom': '4', 'aliTo': '260', 'envFrom': '4', 'envTo': '262', 'acc': '0.99'}]}]
