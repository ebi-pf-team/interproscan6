import argparse
import json


def union_results(pre_calc: list, analysis: list):
    all_results = {}
    results = pre_calc + analysis
    for file_path in results:
        if file_path:
            with open(file_path, 'r') as file:
                data = json.load(file)
            for seq_id, match_info in data.items():
                try:
                    all_results[seq_id].append(match_info)
                except:
                    all_results[seq_id] = match_info
    return all_results


def main():
    parser = argparse.ArgumentParser(
        description="sequences parser"
    )
    parser.add_argument(
        "-pre_calc", "--pre_calc_results", type=str, help="list of pre calc files"
    )
    parser.add_argument(
        "-analysis", "--analysis_results", type=str, help="list of analysis files"
    )
    args = parser.parse_args()

    pre_calc_results = args.pre_calc_results.strip("[]").replace(" ", "").split(",")
    analysis_results = args.analysis_results.strip("[]").replace(" ", "").split(",")
    all_results = union_results(pre_calc_results, analysis_results)
    print(json.dumps(all_results))


if __name__ == "__main__":
    main()
