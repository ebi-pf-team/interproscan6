import json
import sys


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
    with open(hmmer_domtbl, "r") as dtbl_f:
        matches = {}
        version = hmmer_domtbl.split("/")[-1].split("_")[0]
        member_db = hmmer_domtbl.split("/")[-1].split("_")[1].split(".")[0]


        for line in dtbl_f.readlines():
            if line.startswith("#"):
                continue
            if line.startswith("[I6-SITES]"):
                if retrieve_sites:
                    continue
                else:
                    break

            info = line.split()

            if line.startswith("[site]") and retrieve_sites:
                # Retrieve site annotations
                site = get_site(info)
                target_key = str(info[1])
                acc_key = str(info[2].split(".")[0])
                signature = matches[target_key][acc_key]
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
                        if (domain["start"] < site_range[0] < domain["end"]) and (
                                domain["start"] < site_range[1] < domain["end"]):
                            try:
                                if site not in signature["locations"][i]["sites"]:
                                    signature["locations"][i]["sites"].append(site)
                            except KeyError:
                                signature["locations"][i]["sites"] = [site]
            else:
                # Retrieve signature data
                target_key = str(info[0])
                acc_key = str(info[4].split(".")[0])
                signature = get_signature(info, member_db, version)
                location = get_domain(info)
                if target_key not in matches:
                    matches[target_key] = {}
                if acc_key not in matches[target_key]:
                    matches[target_key][acc_key] = signature
                    matches[target_key][acc_key]["locations"] = [location]
                else:
                    if location not in matches[target_key][acc_key]["locations"]:
                        matches[target_key][acc_key]["locations"].append(location)
    return matches


def get_signature(info: list[str], member_db: str, version: str) -> dict[str, str]:
    """
    Retrieve data for the full sequence hit against the model:
        These are the data listed under --- full sequence ---

    :param info: list, line split by blankspace
    :param locations: list of dicts, one dict per hit for the signature
        against the query protein sequence
    """
    signature_info = {
        "accession": info[4].split(".")[0],
        "name": info[3],
        "evalue": float(info[6]),
        "score": float(info[7]),
        "qlen": int(info[5]),
        "bias": float(info[8]),
        "member_db": member_db,
        "version": version,
        "model-ac": info[4]
    }
    return signature_info


def get_domain(info: list[str]) -> dict[str, str]:
    """
    Retrieve the data for the specific domain hit:
        These are the data listed under --- this domain ---
    These data are stored under locations in the output JSON,
    and represent all the places where a InterPro signature matched
    the input query sequence

    :param info: list, line split by blankspace
    """
    domain_info = {
        "start": int(info[17]),  # ali coord from
        "end": int(info[18]),   # ali coord to
        "representative": "",
        "hmmStart": int(info[15]),  # hmm coord from
        "hmmEnd": int(info[16]),  # hmm coord to
        "hmmLength": int(info[5]),  # qlen
        "hmmBounds": "",
        "evalue": float(info[12]),  # Independent e-value
        "score": float(info[13]),  # bit score
        "envelopeStart": int(info[19]),  # env coord from
        "envelopeEnd": int(info[20]),  # env coord to
        "postProcessed": ""
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
