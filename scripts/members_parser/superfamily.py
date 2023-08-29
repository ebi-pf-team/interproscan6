import hashlib
import subprocess
from typing import TextIO

import pyhmmer
from Bio import SeqIO
from pyhmmer.easel import Alphabet, SequenceFile
from pyhmmer.plan7 import HMMFile

project_dir = "/Users/lcf/interproscan6"
hmm_superfamily_path = f"{project_dir}/data/superfamily/1.75/hmmlib_1.75"
all_appl = f"{project_dir}/files_test/test_all_appl.fasta"

bin_ass3_pl_path = f"{project_dir}/bin/superfamily/1.75/ass3_single_threaded.pl"
self_hits_path = f"{project_dir}/data/superfamily/1.75/self_hits.tab"
cla_path = f"{project_dir}/data/superfamily/1.75/dir.cla.scop.txt_1.75"
model_tab_path = f"{project_dir}/data/superfamily/1.75/model.tab"
pdbj95d_path = f"{project_dir}/data/superfamily/1.75/pdbj95d"
hmmer3_hmmsearch_cpu_switch = "--cpu 1"
hmmer3_hmmsearch_switches = "-E 10 -Z 15438"


def get_sequences(fasta_path: str) -> dict:
    sequences = {}
    for sequence in SeqIO.parse(fasta_path, "fasta"):
        sequence_info = []
        sequence_info.append(sequence.name)
        sequence_info.append(sequence.description)
        sequence_info.append(sequence.seq)
        sequence_info.append(hashlib.md5(str(sequence.seq).encode()).hexdigest())
        sequences[sequence.id] = sequence_info
    return sequences


# def _write_full_length_results(
#     fdst: TextIO,
#     name: str,
#     hmmname: str,
#     evalue: float,
#     bitscore: float,
#     dom_evalue: float,
#     dom_bitscore: float,
#     dom_start: int,
#     dom_end: int,
# ):
#     msg = f"{name}\t{hmmname}\t{evalue:.2e}\t{bitscore:.2f}\n"
#     fdst.write(msg)
#
#
# def write_output(allhits, fdst):
#     # if domain_search:
#     #     header = "query\tsearch_homologous\tdom_e-value\tdom_bitscore\tstart\tend\n"
#     #     writer = _write_domain_results
#     # else:
#     header = "query\tsearch_homologous\te-value\tbitscore\n"
#     writer = _write_full_length_results
#
#     fdst.write(header)
#     for hits in allhits:
#         for hit in hits.reported:
#             name = hit.name.decode()
#             hmmname = hits.query_name.decode()
#             evalue = hit.evalue
#             bitscore = hit.score
#             dom_evalue = hit.best_domain.i_evalue
#             dom_bitscore = hit.best_domain.score
#             dom_start = hit.best_domain.env_from
#             dom_end = hit.best_domain.env_to
#             writer(
#                 fdst=fdst,
#                 name=name,
#                 hmmname=hmmname,
#                 evalue=evalue,
#                 bitscore=bitscore,
#                 dom_evalue=dom_evalue,
#                 dom_bitscore=dom_bitscore,
#                 dom_start=dom_start,
#                 dom_end=dom_end,
#             )


def run_hmm(fasta_file: str, hmm_file_path: str, output_path: str):
    alphabet = Alphabet.amino()
    with SequenceFile(fasta_file, digital=True, alphabet=alphabet) as sf:
        sequencesD = sf.read_block()

    with HMMFile(hmm_file_path) as hmm_file:
        all_hits = pyhmmer.hmmscan(sequencesD, hmm_file, E=10, Z=15438, cpus=1)

    for hit in all_hits:
        print(hit.query_name)

    return all_hits


def run_pl_script(hmm_result, output_path):
    command = [
        "perl",
        hmm_superfamily_path,
        hmmer3_hmmsearch_cpu_switch,
        hmmer3_hmmsearch_switches,
        "-s",
        self_hits_path,
        "-r",
        cla_path,
        "-m",
        model_tab_path,
        "-p",
        pdbj95d_path,
        all_appl,
        hmm_result,
        output_path,
    ]

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(e)


if __name__ == "__main__":
    output_pl_path = (
        "/Users/lcf/PycharmProjects/interproscan6/results/output_superfamily_pl.txt"
    )
    output_hmmscan_path = "/Users/lcf/PycharmProjects/interproscan6/results/output_superfamily_hmm_scan.txt"
    sequences = get_sequences(all_appl)
    hmm_result = run_hmm(all_appl, hmm_superfamily_path, output_hmmscan_path)
    # run_pl_script(hmm_result, output_pl_path)
