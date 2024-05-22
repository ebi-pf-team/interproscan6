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
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(ANF\d{5})\s*$")
    if appl.upper() == "NCBIFAM":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+((TIGR|NF)\d+)\.\d+$")
    if appl.upper() == "PANTHER":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(PTHR[^\s]+)\s*\[M=\d+\]$")
    if appl.upper() == "SFLD":
        return re.compile(r"^(Accession:|Query:|Query sequence:)\s+(SFLD[^\s]+)\s*$")


def parse(out_file: str) -> dict:
    version = out_file.split("/")[-1].split("_")[0]
    member_db = out_file.split("/")[-1].split("_")[1].split(".")[0]
    current_sequence = None
    current_domain = None
    domain_match = {}
    hmmer_parser_support = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    appl = out_file.split("_")[1].split(".")[0]
    member_accession = get_accession_regex(appl)
    description = ""

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
                            
                            try:
                                cigar_alignment = cigar_alignment_parser(domain_match[domain_key]["alignment"])
                                domain_match[domain_key]["cigar_alignment"] = encode(cigar_alignment)
                            except KeyError:
                                pass # temp until fixed
                            
                            if "locations" not in sequence_match:
                                sequence_match["locations"] = []
                            sequence_match["locations"].append(domain_match[domain_key])
                        domain_match = {}
                        if current_sequence in hmmer_parser_support:
                            hmmer_parser_support[current_sequence].update({model_id: sequence_match})
                        else:
                            hmmer_parser_support[current_sequence] = {model_id: sequence_match}
                    sequence_match = {}
                    description = ""
                    stage = "LOOKING_FOR_METHOD_ACCESSION"
                else:
                    if stage == 'LOOKING_FOR_METHOD_ACCESSION':
                        if line.startswith(("Accession:", "Query:", "Query sequence:")):
                            stage = 'LOOKING_FOR_SEQUENCE_MATCHES'
                            model_ident_pattern = member_accession.match(line)
                            if model_ident_pattern:
                                model_id = model_ident_pattern.group(2).replace(".orig.30.pir", "")
                        if line.startswith(("Accession:", "Query:", "Query sequence:")):
                            query_name = line.split()[1].replace(".orig.30.pir", "")
                            qlen = line.split("[")[1].split("]")[0].replace("M=", "")
                    elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                        if line.startswith("Description:"):
                            description = line.replace("Description:", "")
                        if line.strip() == "":
                            stage = 'LOOKING_FOR_DOMAIN_SECTION'
                            current_domain = None
                            current_sequence = None
                        else:
                            sequence_match = get_sequence_match(
                                line,
                                model_id,
                                query_name,
                                description,
                                version,
                                member_db,
                                qlen
                            )
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
                                domain_match[domain_number] = get_domain_match(match, member_db, qlen)

    return hmmer_parser_support


def get_domain_match(match: re.Match, member_db: str, qlen: str) -> dict:
    post_processed = "false"
    if member_db == "gene3d" or member_db == "pfam":
        post_processed = "true"
    domain_match = {}
    hmm_bound_pattern = {"[]": "Complete", "[.": "N-terminal complete", ".]": "C-terminal complete", "..": "Incomplete"}
    domain_match["start"] = match.group(9)  # ali coord from
    domain_match["end"] = match.group(10)  # ali coord to
    domain_match["representative"] = ""
    domain_match["hmmStart"] = match.group(6)  # hmm coord from
    domain_match["hmmEnd"] = match.group(7)   # hmm coord to
    domain_match["hmmLength"] = qlen  # qlen
    domain_match["rawHmmBounds"] = match.group(8)
    domain_match["hmmBounds"] = hmm_bound_pattern[match.group(8)]
    domain_match["evalue"] = match.group(5)  # Independent e-value
    domain_match["score"] = match.group(2)   # bit score
    domain_match["envelopeStart"] = match.group(11)    # env coord from
    domain_match["envelopeEnd"] = match.group(12)  # env coord to
    domain_match["postProcessed"] = post_processed

    return domain_match


def get_sequence_match(sequence_line: str, model_id: str, query_name: str, description: str, version: str, member_db: str, qlen: str) -> dict:
    sequence_match = {}
    SEQUENCE_LINE_PATTERN = re.compile("^\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\S+\\s+\\d+\\s+(\\S+).*$")
    match = SEQUENCE_LINE_PATTERN.match(sequence_line)
    if match:
        sequence_match["accession"] = model_id
        sequence_match["name"] = query_name
        sequence_match["description"] = description.strip()
        sequence_match["evalue"] = match.group(1)
        sequence_match["score"] = match.group(2)
        sequence_match["qlen"] = qlen
        sequence_match["bias"] = match.group(3)
        sequence_match["member_db"] = member_db
        sequence_match["version"] = version
        sequence_match["model-ac"] = model_id.split(":")[0].split(".")[0]

    # if member_db.lower() == 'panther':
    #     node_id = info[-1]  # retrieve node id from Panther-TreeGrafter hits
    #     paint_anno_path = mem_db_dir + f'/{info[4].split(":")[0].split(".")[0]}.json'
    #     with open(paint_anno_path, 'r') as fh:
    #         paint_annotations = json.load(fh)
    #         node_data = paint_annotations[node_id]
    #     signature_info['proteinClass'] = node_data[2]
    #     signature_info['graftPoint'] = node_data[3]

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
