"""Run pfsearchV3 for all PROSITE profiles"""
import os
import sys
import subprocess


def get_profile_paths(directory):
    for dirpath, _, filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))


def main():
    if len(sys.argv) < 5:
        sys.exit("Error: expected more than 5 arguments, check your command again")

    models_dir = sys.argv[1]
    fasta_file = sys.argv[2]
    output_file = sys.argv[3]
    bin_cmd = sys.argv[4]  # path to the pfsearch binary or executable command
    binary_switches = sys.argv[5:]

    profiles = get_profile_paths(models_dir)

    for profile in profiles:
        run_cmd = [bin_cmd, profile, fasta_file] + binary_switches
        try:
            output = subprocess.check_output(run_cmd, universal_newlines=True)
            if output.strip():
                with open(output_file, 'a') as out_file:
                    out_file.write(output + '\n')
        except subprocess.CalledProcessError as e:
            sys.exit(f"Error running pfsearchV3 using cmd: {' '.join(run_cmd)}\nError: {e}")


if __name__ == "__main__":
    main()
