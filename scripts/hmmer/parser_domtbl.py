import json
import sys

COMMENT_LINE = "#"


def parse(hmmer_domtbl: str):
    with open(hmmer_domtbl, "r") as dtbl_f:
        matches = {}
        signatures = {}  # accession is the unique key
        locations = []
        member_db = hmmer_domtbl.split("/")[-1].split(".")[0]

        for line in dtbl_f.readlines():
            if line.startswith(COMMENT_LINE):
                continue
            info = line.split()
            target_key = str(info[0])
            acc_key = str(info[4])
            signature = get_signature(info, member_db)
            location = get_domain(info)

            if signatures.get(acc_key) is None:
                signatures[acc_key] = signature
                locations = []
            else:
                locations.append(location)
                signatures[acc_key]["locations"] = locations

            try:
                matches[target_key].append(signatures)
            except KeyError:
                matches[target_key] = [signatures]

    return matches


def get_signature(info, member_db) -> dict:
    signature_info = {
        "accession": info[4],
        "name": info[3],
        "description": "",   # just in hmm.out?
        "e_value": float(info[6]),
        "score": float(info[7]),
        "member_db": member_db
    }
    return signature_info


def get_domain(info) -> dict:
    domain_info = {
        "start": int(info[17]),
        "end": int(info[18]),
        "representative": "",
        "hmmStart": int(info[15]),
        "hmmEnd": int(info[16]),
        "hmmLength": int(info[5]),  # qlen?
        "hmmBounds": "",
        "evalue": float(info[12]),  # iEvalue?
        "score": float((info[13])),
        "envelopeStart": int(info[19]),
        "envelopeEnd": int(info[20]),
        "postProcessed": ""
    }
    return domain_info


def main():
    args = sys.argv[1:]
    parse_result = parse(args[0])
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
