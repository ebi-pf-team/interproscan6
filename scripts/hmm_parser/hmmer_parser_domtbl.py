import argparse
import json

COMMENT_LINE = "#"


def parse(domtbl_file: str):
    sequence_matches = {}
    with open(domtbl_file, "r") as dtbl_f:
        current_seq = None
        acc = []
        domains = []
        for line in dtbl_f.readlines():
            if not line.startswith(COMMENT_LINE):
                info = line.split()
                if info[0] != current_seq:
                    if current_seq:
                        sequence_matches[current_seq] = acc
                    acc = []
                    current_seq = info[0]
                if info[9] == info[10]:
                    domains.append(get_domain(info))
                    acc.append(get_accession(info, domains))
                    domains = []
                else:
                    domains.append(get_domain(info))

            if current_seq:
                sequence_matches[current_seq] = acc
    return sequence_matches


def get_accession(info, domains):
    acc_info = {
        "query_name": info[3],
        "accession": info[4],
        "qlen": int(info[5]),
        "e_value": float(info[6]),
        "score": float(info[7]),
        "bias": float(info[8]),
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
        "accession": info[21],
        "description_of_target": info[22]
    }
    return domain_info


def main():
    parser = argparse.ArgumentParser(
        description="hmmer domtbl parser"
    )

    parser.add_argument(
        "-hmmer_file", "--hmmer_file", type=str, help="dtbl file result of hmmer preproc")
    args = parser.parse_args()

    parse_result = parse(args.hmmer_file)
    print(json.dumps(parse_result, indent=2))


if __name__ == "__main__":
    main()
