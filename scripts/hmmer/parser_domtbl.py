import json
import sys

COMMENT_LINE = "#"


def parse(hmmer_domtbl: str):
    sequence_matches = {}
    with open(hmmer_domtbl, "r") as dtbl_f:
        acc = []
        domains = []
        member_db = hmmer_domtbl.split("/")[-1].split(".")[0].split("_")[1]

        for line in dtbl_f.readlines():
            if line.startswith(COMMENT_LINE):
                continue

            info = line.split()
            seq_id = info[0]

            try:
                acc_info = get_accession(info, domains, member_db)
                domain[] = get_domain(info)
                sequence_matches[seq_id].append()
            if seq_id != current_seq:
                if current_seq:
                    sequence_matches[current_seq] = acc
                acc = []
                current_seq = seq_id
            if info[9] == info[10]:
                domains.append(get_domain(info))
                acc.append(get_accession(info, domains, member_db))
                domains = []
            else:
                domains.append(get_domain(info))

        if current_seq:
            sequence_matches[current_seq] = acc

    return sequence_matches


def get_accession(info, domains, member_db):
    acc_info = {
        "query_name": info[3],
        "accession": info[4].split(".")[0],
        "qlen": int(info[5]),
        "e_value": float(info[6]),
        "score": float(info[7]),
        "bias": float(info[8]),
        "member_db": member_db,
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
    args = sys.argv[1:]
    parse_result = parse(args[0])
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
