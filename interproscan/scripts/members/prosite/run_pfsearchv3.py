"""Run pfsearchV3 for all PROSITE profiles"""
import os
import sys
import subprocess


def get_profile_paths(directory):
    for dirpath, _, filenames in os.walk(directory):
        for f in filenames:
            yield os.path.abspath(os.path.join(dirpath, f))


def main():
    """
    CL input:
    0. Str repr of path to PROSITE profiles model directory
    1. Str repr of path for fasta file
    2. Str repr of path for the output file
    3. Path to the pfsearch binary or executable command
    4. Str repr of path to the IPS6 binary switches
    """
    args = sys.argv[1:]
    models_dir = args[0]
    fasta_file = args[1]
    output_file = args[2]
    bin_cmd = args[3]
    binary_switches = args[4:]

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
