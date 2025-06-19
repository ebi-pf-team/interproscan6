# The superfamily unit tests need the HMM file but it is huge!
# This script strips back the hmm file to a minimal working file for the unit tests

import re

def parse_hmm_file(input_file, output_file):
    model2length = {}
    modelAc = None
    length = None

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            line = line.strip()
            if line.startswith('//'):
                assert modelAc is not None and length is not None
                model2length[modelAc] = length
                outfile.write(f"NAME {modelAc}\n")
                outfile.write(f"LENG {length}\n")
                outfile.write("//\n")
                modelAc = length = None
            elif line.startswith('N') and not modelAc:
                match = re.match(r'^NAME\s+(.+)$', line)
                if match:
                    modelAc = match.group(1)
            elif line.startswith('L') and not length:
                match = re.match(r'^LENG\s+([0-9]+)$', line)
                if match:
                    length = int(match.group(1))

parse_hmm_file('tests/data/databases/superfamily/1.75/hmmlib_1.75', 'tests/data/databases/superfamily/1.75/hmmlib_1.75_minimalist')
