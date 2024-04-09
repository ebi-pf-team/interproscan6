import json
import sys


def aggregate_results(result_files: list):
    all_results = {}
    for file_path in result_files:
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
    args = sys.argv[1:]
    result_files = args[0].strip('[]').replace(" ", "").split(',')
    all_results = aggregate_results(result_files)
    print(json.dumps(all_results, indent=4))


if __name__ == "__main__":
    main()
