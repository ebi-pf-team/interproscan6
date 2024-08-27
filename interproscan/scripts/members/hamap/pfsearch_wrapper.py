"""Script to parse the hmmer3 tbl file and run the pfsearch binary

It does not filter the matches (reports only normal (0) amd strong (1) matches) 
We capture the flag but do not need to check it.
But there are some PROSITE profiles that are considered are ambiguous,
generating too many FPs, these profiles are black listed and all
hits are ignored.
"""


import subprocess
import sys
import os.path
import os
import re
import time


def get_hamap_profile(hmmer_tbl_path: str, model_dir: str) -> dict:
    """Get the hamap profiles from the hmmer3 .tbl file"""
    profiles = {}
    pattern = re.compile(r'^(\S+)\s+(\S+)\s+(\S+)\s+(.*)')
    temp_err_file = hmmer_tbl_path + '-filtered'

    with open(hmmer_tbl_path, "r") as profile_list:
        for line in profile_list:
            line = line.strip()

            if not line.startswith('#'):
                m = pattern.match(line)

                if m:
                    seq_id, _, profile, _ = m.groups()
                    profile_path = model_dir + '/' + profile + ".prf"

                    if profile in profiles:
                        profiles[profile].append(seq_id)
                    else:
                        profiles[profile] = [profile_path, seq_id]

                else:
                    sys.stderr.write(f'Something wrong in hmmer.tbl file.\nLine: {line}')
                    raise ValueError

    key_count = len(profiles)
    output = f'profile hits: {key_count}\n'
    with open(temp_err_file, 'a') as out_file:
        out_file.write(output)

    return profiles


def get_sequences(fasta_file: str) -> dict:
    """Build a dict of keyed by seq id, valued by protein sequences"""
    fasta_dict = {}
    with open(fasta_file, 'r') as fasta:
        for line in fasta:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                seq_id = line.lstrip('>').split(" ")[0]
                fasta_dict[seq_id] = ''
            else:
                fasta_dict[seq_id] += line + '\n'
    return fasta_dict


def run_pfsearch_binary(
    hmmer_tbl_path: str,
    output_file: str,
    profiles,
    seqs_dict: dict,
    input_fasta_file: str,
    fasta_filtered_outfile: str,
    pfsearch_flags: list[str],
    pfsearch_bin
):
    """pfsearch flags:
    -f: input sequence file is in FASTA format.
    --output-method <uint>  [-o] : printing output method
    == 7 tsv output (single line tab delimited)
    """
    def _append_to_file(filename, output):
        with open(filename, 'a') as out_file:
            out_file.write(output)

    def _get_sequences_for_profile(key_list, seqs_dict):
        sequences = ''
        for key in key_list:
            if key in seqs_dict:
                value = seqs_dict[key]
                sequences += '>' + key + '\n' + value
                if not value.endswith('\n'):
                    sequences += '\n'
            else:
                print("key not found : " + key)
        print("KEY_LIST::", key_list)
        print("SEQS DICT::", seqs_dict)
        return sequences

    def _write_to_file(filename, output):
        with open(filename, 'w') as seq_file:
            seq_file.write(output)

    count = 0
    temp_file_list = []
    stats_filename = input_fasta_file + ".stats"
    temp_dir = input_fasta_file + '-tmp'
    os.makedirs(temp_dir)
    temp_err_file = hmmer_tbl_path + '-filtered'

    _append_to_file(temp_err_file, 'profile hits in run_pfsearch: ' + str(len(list(profiles.keys()))) + '\n')

    _append_to_file(temp_err_file, 'sequence count in run_pfsearch: ' + str(len(list(seqs_dict.keys()))) + '\n')
    get_seq_time = 0
    keys = list(profiles.keys())
    _write_to_file(temp_err_file + '-keys', 'profile keys in run_pfsearch: ' + str(len(keys)) + '\n')

    for prf in profiles:
        print("CRNT PRF::", prf)
        prf_seqs = ' '.join(profiles[prf])
        prf_seqs_count = len(profiles[prf][1:])
        prf_out = 'processing profile #:' + str(count) + ' ' +  prf + ' seqs:' + str(prf_seqs_count) + ' - ' + prf_seqs + '  \n'
        _append_to_file(temp_err_file, prf_out)

        if prf == 'MF_00005':   # black listed?
            mf_output = ' '.join(profiles[prf])

        sequence_ids = profiles[prf][1:]
        get_seq_start_time = time.time()
        input_fasta_sequences = _get_sequences_for_profile(sequence_ids, seqs_dict)
        temp_file = temp_dir + '/' + prf

        _write_to_file(temp_file, input_fasta_sequences)
        _append_to_file(fasta_filtered_outfile, input_fasta_sequences)
        get_seq_end_time = time.time()
        time_to_get_seqences = get_seq_end_time - get_seq_start_time

        get_seq_time += time_to_get_seqences
        temp_file_list.append(temp_file)

        comd_to_run = [pfsearch_bin, profiles[prf][0], temp_file] + pfsearch_flags
        count += 1

        if not os.path.isfile(profiles[prf][0]):
            # profile not available
            continue

        output = subprocess.check_output(comd_to_run, universal_newlines=True)
        if output.strip():
            with open(output_file, 'a') as out_file:
                out_file.write(output + '\n')

    _append_to_file(temp_err_file, 'completed running thru ' + str(count) + ' profiles \n')
    with open(stats_filename, 'w') as stats_file:
        stats_file.write('Total time to get and write ' + str(len(temp_file_list)) + ' seq files :' + str(get_seq_time * 1000) + " ms \n")
    return count


def main():
    arg_list = sys.argv
    hmmer_tbl_path = sys.argv[1]
    fasta_file = sys.argv[2]
    fasta_filtered_outfile = sys.argv[3]
    output_file = sys.argv[4]
    model_dir = sys.argv[5]
    pfsearch_bin = sys.argv[6]
    pfsearch_flags = sys.argv[7:]  # something like ['-f', '-o', '7']

    # create the output file in case we don't have any matches
    open(output_file, 'a').close()
    open(fasta_filtered_outfile, 'a').close()

    # get the profiles
    profiles = get_hamap_profile(hmmer_tbl_path, model_dir)

    if profiles:
        # get the protein sequences
        seqs_dict = get_sequences(fasta_file)

        # run the pfsearch binary
        pfsearch_cmd_run_count = run_pfsearch_binary(
            hmmer_tbl_path,
            output_file,
            profiles,
            seqs_dict,
            fasta_file,
            fasta_filtered_outfile,
            pfsearch_flags,
            pfsearch_bin
        )
        sys.stderr.write(f'prfs: {str(pfsearch_cmd_run_count)}')


if __name__ == "__main__":
    if len(sys.argv) < 6:
        sys.exit("Error: expected more than 6 arguments, check your command again")

    main()
