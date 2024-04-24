import json
import sys

COMMENT_LINE = "#"


def parse(hmmer_domtbl: str) -> dict:
    with open(hmmer_domtbl, "r") as dtbl_f:
        matches = {}
        member_db = hmmer_domtbl.split("/")[-1].split(".")[0]

        for line in dtbl_f.readlines():
            if line.startswith(COMMENT_LINE):
                continue
            info = line.split()
            target_key = str(info[0])
            acc_key = str(info[4])
            signature = get_signature(info, member_db)
            location = get_domain(info)

            if target_key not in matches:
                matches[target_key] = {}

            if acc_key not in matches[target_key]:
                matches[target_key][acc_key] = signature
                matches[target_key][acc_key]["locations"] = [location]
            else:
                matches[target_key][acc_key]["locations"].append(location)

    return matches


def get_signature(info, member_db) -> dict:
    signature_info = {
        "accession": info[4],
        "name": info[3],
        "description": "",   # just in hmm.out?
        "e_value": float(info[6]),
        "score": float(info[7]),
        "bias": float(info[8]),
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
        "bias": info[14],
        "postProcessed": ""
    }
    return domain_info


def main():
    args = sys.argv[1:]
    parse_result = parse(args[0])
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
