import json
import re
import sys
"""Parse the output from HMMER3 hmmscan into the internal IPS6 JSON structure.

Note the PIRSF uses the envelope start and stop for its start and end locations,
respectively, instead of 'ali from' and 'ali to' like other member db."""


ACC_LINE = re.compile(r"^Query:\s+(\S+)\s+\[[L|M]=(\d+)\]$")
NO_HITS_LINE = "[No hits detected that satisfy reporting thresholds]"
MODEL_HIT_LINE = re.compile(r"^([\w+\.\-]+)\s+(\d+|[\d\.]+)\s+(\d+|[\d\.]+)\s+([\w+\.\-]+)\s+(\d+|[\d\.]+)\s+(\d+|[\d\.]+)\s+(\d+|[\d\.]+)\s+(\d+)\s+(PIRSF\d+)$")
DOMAIN_HIT_LINE = re.compile(r"^(\d+)\s+(!|\?)\s+([\d\.]+)\s+([\d\.]+)\s+([\d\w\.-]+)\s+([\d\w\.-]+)\s+(\d+)\s+(\d+)\s+([\.\[\]]{2})\s+(\d+)\s+(\d+)\s+([\.\[\]]{2})\s+(\d+)\s+(\d+)\s+([\.\[\]]{2})\s+([\d\.]+)$")
HMM_BOUND_PATTERN = {
    "[]": "COMPLETE",
    "[.": "N_TERMINAL_COMPLETE",
    ".]": "C_TERMINAL_COMPLETE",
    "..": "INCOMPLETE"
}


class QueryProtein:
    """Represents an input protein sequence and its associated hits"""
    def __init__(self):
        self.sequence_id = None
        self.qlen = None
        self.signatures = {}  # sig acc: model hit

    def get_seq_id(self, re_acc_line_match):
        self.sequence_id = re_acc_line_match.group(1)
        self.qlen = re_acc_line_match.group(2)

    def get_model_data(self, re_match):
        """Get data on a signature/model that matched.
        Where model_pattern are the groups from the line matching SIG_DATA_LINE"""
        model = ModelHit()
        model.model_id = re_match.group(9)
        model.evalue = re_match.group(1)
        model.score = re_match.group(2)
        model.bias = re_match.group(3)
        model.best_domain_evalue = re_match.group(4)
        model.best_domain_score = re_match.group(5)
        model.best_domain_bias = re_match.group(6)
        model.num_domain_exp = re_match.group(7)
        model.num_domains = re_match.group(8)
        self.signatures[model.model_id] = model

    def get_domin_data(self, model_id, re_match):
        """Get stats for a domain hit that matched.
        Where domain_pattern are the groups from the line matching DOMAIN_DATA_LINE"""
        domain = DomainHit()
        domain.model_id = model_id
        domain.domain_num = re_match.group(1)
        domain.score = re_match.group(3)
        domain.bias = re_match.group(4)
        domain.c_evalue = re_match.group(5)
        domain.i_evalue = re_match.group(6)
        domain.hmm_from = re_match.group(7)
        domain.hmm_to = re_match.group(8)
        domain.ali_from = re_match.group(10)
        domain.ali_to = re_match.group(11)
        domain.env_from = re_match.group(13)
        domain.env_to = re_match.group(14)
        # In the hmmsearch parser for the other member dbs, we use the hmm_bounds
        # from between the hmm coord and ali coord
        # But for PIRSF we use the bounds associated with the envelope
        # This was something decided by the PIRSF devs
        domain.hmm_raw_bounds = re_match.group(15)
        domain.hmm_bounds = HMM_BOUND_PATTERN[domain.hmm_raw_bounds]
        domain.acc = re_match.group(16)
        self.signatures[domain.model_id].domains[domain.domain_num] = domain


class ModelHit:
    """Represent a signature"""
    def __init__(self):
        self.model_id = None
        self.evalue = None
        self.score = None
        self.bias = None
        self.best_domain_evalue = None
        self.best_domain_score = None
        self.best_domain_bias = None
        self.num_domain_exp = None
        self.num_domains = None
        self.domains = {}  # {domain num: DomainHit}


class DomainHit:
    """Present a domain hit for a signature against a query protein sequence"""
    def __init__(self):
        self.model_id = None
        self.domain_num = None
        self.score = None
        self.bias = None
        self.c_evalue = None
        self.i_evalue = None
        # Below: the start and end points with respect to the consensus coordinates of the model
        self.hmm_from = None
        self.hmm_to = None
        # Below: the start and end points of the alignment on the target sequence.
        self.ali_from = None
        self.ali_to = None
        self.env_from = None
        self.env_to = None
        # Hmm bounds: alignment ended internally in the hmm, or ran all the way to the end
        self.hmm_raw_bounds = None
        self.hmm_bounds = None
        self.acc = None


def add_match(
    matches: dict,
    protein_with_hit: QueryProtein,
    member_db: str,
    version: str,
) -> dict:
    """Store data for protein with hits against SMART HMM profiles.

    Keyed by query protein ID, into dict:
    Keyed by signature accession, into dict:
        accession, name, descriptiion, etc. locations: [],"""
    if protein_with_hit.sequence_id not in matches:
        matches[protein_with_hit.sequence_id] = {}

    for model_id, model_obj in protein_with_hit.signatures.items():  # dict {model id: ModelHit()}
        if model_id not in matches[protein_with_hit.sequence_id]:
            matches[protein_with_hit.sequence_id][model_id] = {
                "accession": model_id,
                "name": "",
                "description": "",
                "qlen": int(protein_with_hit.qlen),
                "evalue": model_obj.evalue,
                "score": model_obj.score,
                "member_db": member_db,
                "version": version,
                "model-ac": model_id,
                "locations": []
            }

        for domain_num, domain_obj in model_obj.domains.items():
            # the PIRSF perl script uses the envelope start/end
            matches[protein_with_hit.sequence_id][model_id]["locations"].append(
                {
                    "start": int(domain_obj.ali_from),
                    "end": int(domain_obj.ali_to),
                    "representative": "",
                    "hmmStart": int(domain_obj.hmm_from),
                    "hmmEnd": int(domain_obj.hmm_to),
                    "hmmLength": int(protein_with_hit.qlen),
                    "rawHmmBounds": domain_obj.hmm_raw_bounds,
                    "hmmBounds": domain_obj.hmm_bounds,
                    "evalue": domain_obj.i_evalue,  # keep as str because can be Xe-Y
                    "score": domain_obj.score,
                    "envelopeStart": int(domain_obj.env_from),
                    "envelopeEnd": int(domain_obj.env_to),
                    "location-fragments": [
                        {
                            "start": int(domain_obj.env_from),
                            "end": int(domain_obj.env_to),
                            "dc-status": "CONTINUOUS",
                        }
                    ]
                }
            )

    return matches


def parse(hmmer_out_path: str):
    """Parse the output from HMMER3 hmmscan into a dict

    :param hmmer_out_path: str repr of path to hmmscan out file
    """
    matches = {}
    version = hmmer_out_path.split("/")[-1].split("._.")[0]
    member_db = hmmer_out_path.split("._.")[1]
    query_protein = QueryProtein()
    current_model = None
    stage = "GET_QUERY_PROTEIN"

    with open(hmmer_out_path, 'r') as reader:
        for line in reader:

            if line.startswith("//"):
                if query_protein.signatures:
                    matches = add_match(matches, query_protein, member_db, version)
                # start a new protein instance
                query_protein = QueryProtein()
                stage = 'LOOKING_FOR_METHOD_ACCESSION'

            elif line.strip() == NO_HITS_LINE:
                stage = "GET_QUERY_PROTEIN"

            elif line.strip().startswith("Query:"):
                query_protein.get_seq_id(ACC_LINE.match(line.strip()))
                stage = "LOOKING_FOR_METHOD_ACCESSION"

            elif stage == 'LOOKING_FOR_METHOD_ACCESSION':
                model_data_line = MODEL_HIT_LINE.match(line.strip())
                if model_data_line:
                    query_protein.get_model_data(model_data_line)
                elif not line.strip():
                    stage = "LOOKING_FOR_DOMAIN_SECTION"

            elif line.startswith(">>"):
                stage = "LOOKING_FOR_DOMAIN_SECTION"
                current_model = line.split()[1]

            elif stage == "LOOKING_FOR_DOMAIN_SECTION":
                if line.strip():
                    stage = "LOOKING_FOR_DOMAIN_SECTION"

                    domain_line_data = DOMAIN_HIT_LINE.match(line.strip())
                    if domain_line_data:
                        query_protein.get_domin_data(current_model, domain_line_data)

    return matches


def main():

    parse_result = parse(sys.argv[1])
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
