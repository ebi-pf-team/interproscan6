import json
import re
import sys
"""Parse the output from HMMER2 hmmpfam for SMART

For a description of the Hmmpfam output file structure see the HMMER-2.3.2 UserGuide
under section 'parsing the domain structure of a sequence with hmmpfam' on 
page 26."""


ALN_DOMAIN_LINE = re.compile(r'\S+:\s+domain \d+ of \d+,')
SIGNATURE_DATA_LINE = re.compile(r"^(\w+?)\s+([\w\s\-\|,\[\]\(\)\.]+?)\s+(\d[e0-9\-\.]+?)\s+([e0-9_\-.]+?)\s+(\d+?)$")
DOMAIN_DATA_LINE = re.compile(r"^(\w+)\s+(\d+/\d+)\s+(\d+)\s+(\d+)\s+([\.\[\]]{2})\s+(\d+)\s+(\d+)\s+([\.\[\]]{2})\s+(\d+|\d+\.\d+)\s+(.+)")
HMM_BOUND_PATTERN = {
    "[]": "COMPLETE",
    "[.": "N_TERMINAL_COMPLETE",
    ".]": "C_TERMINAL_COMPLETE",
    "..": "INCOMPLETE"
}


class QueryProtein:
    """Represents a input protein sequence and associated hits"""
    def __init__(self):
        self.sequence_id = None  # query seq id
        self.signatures = {}  # {sig acc/model id: ModelHit}
        self.has_hits = False

    def get_seq_id(self, value):
        """Get query protein sequence ID. Where value is a line"""
        self.sequence_id = " ".join(value.split()[2:])

    def get_model_data(self, model_pattern):
        """Get data on a signature/model that matched.
        Where model_pattern are the groups from the line matching SIG_DATA_LINE"""
        self.has_hits = True
        model = ModelHit()
        model.model_id = model_pattern.group(1)
        model.score = model_pattern.group(3)
        model.evalue = model_pattern.group(4)
        model.num_domains = int(model_pattern.group(5))
        self.signatures[model.model_id] = model

    def get_domain_data(self, domain_pattern):
        """Get stats for a domain hit that matched.
        Where domain_pattern are the groups from the line matching DOMAIN_DATA_LINE"""
        domain = DomainHit()
        domain.model_id = domain_pattern.group(1)
        domain.domain_num = domain_pattern.group(2).split("/")[0]  # e.g. 1/2
        domain.seq_from = domain_pattern.group(3)
        domain.seq_to = domain_pattern.group(4)
        domain.hmm_from = domain_pattern.group(6)
        domain.hmm_to = domain_pattern.group(7)
        domain.hmm_raw_bounds = domain_pattern.group(8)
        domain.hmm_bounds = HMM_BOUND_PATTERN[domain_pattern.group(8)]
        domain.score = domain_pattern.group(9)
        domain.evalue = domain_pattern.group(10)
        self.signatures[domain.model_id].domains[domain.domain_num] = domain


class ModelHit:
    """Represent a SMART signature"""
    def __init__(self):
        self.model_id = None
        self.score = None
        self.evalue = None
        self.num_domains = None
        self.domains = {}  # {domain num: DomainHit}


class DomainHit:
    """Present a domain hit for a signature against a query protein sequence"""
    def __init__(self):
        self.model_id = None
        self.domain_num = None
        # Below: the start and end points of the alignment on the target sequence.
        self.seq_from = None
        self.seq_to = None
        # Below: the start and end points with respect to the consensus coordinates of the model
        self.hmm_from = None
        self.hmm_to = None
        # Hmm bounds: alignment ended internally in the hmm, or ran all the way to the end
        self.hmm_raw_bounds = None
        self.hmm_bounds = None
        self.score = None
        self.evalue = None


def add_match(
    matches: dict,
    protein_with_hit: QueryProtein,
    member_db: str,
) -> dict:
    """Store data for protein with hits against SMART HMM profiles.

    Keyed by query protein ID, into dict:
    Keyed by signature accession, into dict:
        accession, evalue, score, etc. locations: [],"""
    if protein_with_hit.sequence_id not in matches:
        matches[protein_with_hit.sequence_id] = {}

    for model_id, model_obj in protein_with_hit.signatures.items():  # dict {model id: ModelHit()}
        if model_id not in matches[protein_with_hit.sequence_id]:
            matches[protein_with_hit.sequence_id][model_id] = {
                "accession": model_id,
                "evalue": model_obj.evalue,
                "score": model_obj.score,
                "member_db": member_db,
                "model-ac": model_id,
                "locations": []
            }

        for domain_num, domain_obj in model_obj.domains.items():
            matches[protein_with_hit.sequence_id][model_id]["locations"].append(
                {
                    "start": int(domain_obj.seq_from),
                    "end": int(domain_obj.seq_to),
                    "hmmStart": int(domain_obj.hmm_from),
                    "hmmEnd": int(domain_obj.hmm_to),
                    "hmmLength": int(domain_obj.hmm_to) + 1 - int(domain_obj.hmm_from),
                    "rawHmmBounds": domain_obj.hmm_raw_bounds,
                    "hmmBounds": domain_obj.hmm_bounds,
                    "evalue": domain_obj.evalue,  # keep as str because can be Xe-Y
                    "score": domain_obj.score,  # keep as str because can be Xe-Y
                    "postProcessed": "false",
                    "location-fragments": [
                        {
                            "start": int(domain_obj.seq_from),
                            "end": int(domain_obj.seq_to),
                            "dc-status": "CONTINUOUS",
                        }
                    ]
                }
            )

    return matches


def parse_hmmpfam_out(out_file: str, member_db: str) -> dict:
    """Coordinate parsing the HMMER2 Hmmpfam output file"""
    matches = {}
    stage = 'LOOKING_FOR_METHOD_ACCESSION'
    protein_with_hit = QueryProtein()

    with open(out_file, "r") as fh:
        for line in fh.readlines():
            if line.startswith("//"):
                # add new domain to matches
                if protein_with_hit.signatures:
                    matches = add_match(matches, protein_with_hit, member_db)

                # start a new protein instance
                protein_with_hit = QueryProtein()
                stage = 'LOOKING_FOR_METHOD_ACCESSION'

            elif stage == 'LOOKING_FOR_METHOD_ACCESSION':
                if line.startswith("Query sequence:"):
                    protein_with_hit.get_seq_id(line)
                    # go onto next step: getting model matches
                    stage = 'LOOKING_FOR_SEQUENCE_MATCHES'

            elif stage == 'LOOKING_FOR_SEQUENCE_MATCHES':
                if line.strip() == '[no hits above thresholds]':
                    # not matches for this protein so go wait till the next seq
                    stage = 'LOOKING_FOR_METHOD_ACCESSION'

                elif line.startswith('Parsed for domains:'):
                    # onto the next step: getting domain hit data
                    stage = 'LOOKING_FOR_DOMAIN_SECTION'

                else:
                    # see if line contains data for signature/model
                    model_line_pattern = SIGNATURE_DATA_LINE.match(line)
                    if model_line_pattern:
                        protein_with_hit.get_model_data(model_line_pattern)

            elif stage == 'LOOKING_FOR_DOMAIN_SECTION':
                if line.startswith('Alignments of top-scoring domains:'):
                    stage = 'LOOKING_FOR_METHOD_ACCESSION'

                else:
                    domain_line_pattern = DOMAIN_DATA_LINE.match(line)
                    if domain_line_pattern:
                        protein_with_hit.get_domain_data(domain_line_pattern)

    return matches


def main():
    """
    :args 0: str repr of path to hmmer file to be parsed
    :args 1: member database
    :args 2: str repr of path to output file
    """
    args = sys.argv[1:]
    parse_result = parse_hmmpfam_out(args[0], args[1])
    with open(args[2], "w") as fh:
        json.dump(parse_result, fh)


if __name__ == "__main__":
    main()
