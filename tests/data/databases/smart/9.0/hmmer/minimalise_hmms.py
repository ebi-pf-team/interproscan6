# We need the HMM files for the unit tests put they can be big
# So this script strips them back to the minimum data we need
# for the unit tests to work

import os

def minimalise_hmm_files(input_dir, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for filename in os.listdir(input_dir):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        if not os.path.isfile(input_path):
            continue

        acc = None
        leng = None

        with open(input_path, 'r') as infile:
            lines = infile.readlines()

        for line in lines:
            if line.startswith("ACC "):
                acc = line
            elif line.startswith("LENG "):
                leng = line
            elif line.startswith("//"):
                break  # Stop after the first model

        if acc and leng:
            with open(output_path, 'w') as outfile:
                outfile.write(acc)
                outfile.write(leng)
                outfile.write("//\n")

minimalise_hmm_files("tests/data/databases/smart/9.0/hmmer", "tests/data/databases/smart/9.0/hmmer")
