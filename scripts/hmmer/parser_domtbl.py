import json
import sys

COMMENT_LINE = "#"


def parse(hmmer_domtbl: str, retrieve_sites: bool):
    """Parse hmmer output into a JSON object.

    :param hmmer_domtbl: str repr of path to hmmer.dtbl file
    :param retrieve_sites: bool, if to retrieve site annotations
        [only true for SFLD and CDD.]
        [Site annotation is added to the hmmer output files by
        in house post-processing scripts.]
    """
    sequence_matches = {}
    with open(hmmer_domtbl, "r") as dtbl_f:
        current_seq = None
        models = []
        domains = []
        for line in dtbl_f.readlines():
            if not line.startswith(COMMENT_LINE):
                info = line.split()

                if line.startswith("[I6-SITES]"):
                    # add last domain hit to be processed
                    sequence_matches[current_seq] = models

                    if retrieve_sites:
                        continue
                    else:
                        break

                if line.startswith("[site]") and retrieve_sites:
                    current_seq = info[1]
                    for _site in get_site(info):
                        for i, model in enumerate(sequence_matches[current_seq]):
                            if model["model_ac"] == _site["model_ac"]:
                                try:
                                    sequence_matches[current_seq][i]["sites"].append(_site)
                                except KeyError:
                                    sequence_matches[current_seq][i]["sites"] = [_site]
                    continue

                if info[0] != current_seq:
                    if current_seq:
                        sequence_matches[current_seq] = models
                    models = []
                    current_seq = info[0]
                if info[9] == info[10]:  # domain number == num of domains
                    domains.append(get_domain_data(info))
                    models.append(get_full_seq_data(info, domains))
                    domains = []
                else:
                    domains.append(get_domain_data(info))

            if current_seq and retrieve_sites is False:
                sequence_matches[current_seq] = models

    return sequence_matches


def get_full_seq_data(info: list[str], domains: dict[str, str]) -> dict[str, str]:
    """
    Retrieve data for the full sequence hit against the model:
        These are the data listed under --- full sequence ---

    :param info: list, line split by blankspace
    :param domains: list of dicts, one dict per domain hit
    """
    model_info = {
        "query_name": info[3],
        "model_ac": info[4],
        "qlen": int(info[5]),
        "e_value": float(info[6]),
        "score": float(info[7]),
        "bias": float(info[8]),
        "domains": domains,
    }
    return model_info


def get_domain_data(info: list[str]) -> dict[str, str]:
    """
    Retrieve the data for the specific domain hit:
        These are the data listed under --- this domain ---

    :param info: list, line split by blankspace
    """
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
        "modelsession": info[21],
        "description_of_target": info[22]
    }
    return domain_info


def get_site(info: list[str]):
    """
    Retrieve data for site hits. SFLD and CDD only.

    :param info: list, line split by blankspace
    """
    for residue in info[3].split(","):
        yield {
            "model_ac": info[2],
            "residue": residue,
            "description": " ".join(info[4:])
        }


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
