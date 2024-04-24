import json
import sys

COMMENT_LINE = "#"


def parse(hmmer_domtbl: str, retrieve_sites: bool):
    """Parse hmmer output into a JSON object.

    In the resulting dict, each query protein is represented by its 
    query protein ID, extracted directly from the input FASTA file.

    The dict is keyed by protein IDs, and valued by lists of matches:
    one match (represented as a dict) per InterPro signature match
    
    :param hmmer_domtbl: str repr of path to hmmer.dtbl file
    :param retrieve_sites: bool, if to retrieve site annotations
        [only true for SFLD and CDD.]
        [Site annotation is added to the hmmer output files by
        in house post-processing scripts.]
    """
    sequence_matches = {}
    with open(hmmer_domtbl, "r") as dtbl_f:
        current_seq = None
        signatures = []
        domains = []
        for line in dtbl_f.readlines():
            if not line.startswith(COMMENT_LINE):
                info = line.split()

                if line.startswith("[I6-SITES]"):
                    # add last domain hit to be processed
                    sequence_matches[current_seq] = signatures

                    if retrieve_sites:
                        continue
                    else:
                        break

                if line.startswith("[site]") and retrieve_sites:
                    current_seq = info[1]
                    site = get_site(info)
                    sig_acc = info[2]
                    # identify the corresponding signature and domain ('location') hit

                    for signature in sequence_matches[current_seq]:
                        if signature["signature_acc"] == sig_acc:
                            if len(signature["locations"]) == 1:
                                try:
                                    if site not in signature["locations"][0]["sites"]:
                                        signature["locations"][0]["sites"].append(site)
                                except KeyError:
                                    signature["locations"][0]["sites"] = [site]
                            else:
                                site_locations = []
                                for _ in site["siteLocations"]:
                                    site_locations.extend([int(_["start"]), int(_["end"])])
                                site_range = (min(site_locations), max(site_locations))
                                for i, domain in enumerate(signature["locations"]):
                                    if (domain["start"] < site_range[0] < domain["end"]) and (domain["start"] < site_range[1] < domain["end"]):
                                        try:
                                            if site not in signature["locations"][i]["sites"]:
                                                signature["locations"][i]["sites"].append(site)
                                        except KeyError:
                                            signature["locations"][i]["sites"] = [site]
                    continue

                if info[0] != current_seq:
                    if current_seq:
                        sequence_matches[current_seq] = signatures
                    signatures = []
                    current_seq = info[0]

                if info[9] == info[10]:
                    # domain number == num of domains, therefore parsed all domain hits for this protein
                    domain = get_domain_hit_data(info)
                    if domain not in domains:
                        domains.append(domain)
                    signatures.append(get_signature_data(info, domains))
                    domains = []

                else:
                    domain = get_domain_hit_data(info)
                    if domain not in domains:
                        domains.append(domain)

            if current_seq and retrieve_sites is False:
                sequence_matches[current_seq] = signatures

    return sequence_matches


def get_signature_data(info: list[str], locations: dict[str, str]) -> dict[str, str]:
    """
    Retrieve data for the full sequence hit against the model:
        These are the data listed under --- full sequence ---

    :param info: list, line split by blankspace
    :param locations: list of dicts, one dict per hit for the signature
        against the query protein sequence
    """
    signature_info = {
        "signature_acc": info[4],
        "signature_name": info[3],
        "evalue": float(info[6]),
        "score": float(info[7]),
        "qlen": int(info[5]),
        "bias": float(info[8]),
        "locations": locations,
    }
    return signature_info


def get_domain_hit_data(info: list[str]) -> dict[str, str]:
    """
    Retrieve the data for the specific domain hit:
        These are the data listed under --- this domain ---
    These data are stored under locations in the output JSON, 
    and represent all the places where a InterPro signature matched
    the input query sequence

    :param info: list, line split by blankspace
    """
    domain_info = {
        "start": info[17],  # ali coord from
        "end": info[18],   # ali coord to
        "hmmStart": info[15],  # hmm coord from
        "hmmEnd": info[16],  # hmm coord to
        "envelopeStart": info[19],  # env coord from
        "envelopeEnd": info[20],  # env coord to
        "cEvalue": info[11],  # Conditional e-value
        "iEvalue": info[12],  # Independent e-value
        "score": info[13],  # bit score
        "bias": info[14],
        "signature_accession": info[21],
        "description_of_target": info[22]
    }
    return domain_info


def get_site(info: list[str]) -> dict:
    """
    Retrieve data for site hits. SFLD and CDD only.

    :param info: list, line split by blankspace
    """
    _site = {
        "description": " ".join(info[4:]),
        "numLocations": len(info[3].split(",")),
        "siteLocations": [],
    }

    for residue in info[3].split(","):
        if residue.find("-") != -1:
            res_start = residue.split("-")[0][1:]
            res_end = residue.split("-")[1]
        else:
            res_start = res_end = residue[1:]
        _site["siteLocations"].append({
            "start": res_start,
            "end": res_end,
            "residue": residue[0]
        })

    return _site


def main():
    """
    :args 0: str repr of path to hmmer file to be parsed
    :args 1: bool whether to retrieve site annotations
        from file [true for SFLD and CDD]
    """
    args = sys.argv[1:]
    parse_result = parse(args[0], True if args[1] == "true" else False)
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
