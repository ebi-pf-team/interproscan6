from Bio import SeqIO
from Bio.Blast.Applications import NcbirpsblastCommandline
from Bio.Blast import NCBIXML
import subprocess

bin_path = "/opt/interproscan6/bin"
project_dir = "/opt/interproscan6"

all_appl = f"{project_dir}/files_test/test_all_appl.fasta"
single_protein = f"{project_dir}/files_test/test_single_protein.fasta"
fasta_path = f"{project_dir}/results"

cdd_signature_library_release = "3.20"
cdd_tbl_path = f"{project_dir}/data/cdd/3.20/data/cddid.tbl"
cdd_library_path = f"{project_dir}/data/cdd/3.20/db/Cdd_NCBI"
cdd_data_path = f"{project_dir}/data/cdd/3.20/data"
binary_rpsblast_path = f"{bin_path}/cdd/rpsblast"
binary_rpsbproc_path = f"{bin_path}/cdd/rpsbproc"
rpsblast_switches_cdd = "-evalue 0.01 -seg no -outfmt 11"  # output json 11, xml 5
rpsbproc_switches_cdd = "-m std"
output_blast_path = f"{project_dir}/results/cdd_blast_output/rpsblast.blast.raw.out"


def get_sequences(fasta_path: str):
    sequences = []
    for record in SeqIO.parse(fasta_path, "fasta"):
        sequences.append(record)


def run_rpsblast():
    cmd = NcbirpsblastCommandline(
        cmd=binary_rpsblast_path,
        query=all_appl,
        db=cdd_library_path,
        evalue=0.01,
        seg="no",
        outfmt=11,
        out=output_blast_path
    )
    cmd()


def run_rpsbproc():
    output_file = f"{project_dir}/results/cdd_blast_output/output_rpsbproc.txt"

    command = [
        binary_rpsbproc_path,
        "-i", output_blast_path,
        "-d", cdd_data_path, "--data-mode", "std"]

    with open(output_file, "w") as out:
        try:
            subprocess.run(command, check=True, stdout=out)
        except subprocess.CalledProcessError as e:
            print(e)


def parse_blast():
    result_handle = open("./results/output_rpsblast.xml")
    blast_records = NCBIXML.parse(result_handle)
    for blast_record in blast_records:
        for alignment in blast_record.alignments:
            print("Hit_id:", alignment.hit_id)
            print("Hit_def:", alignment.hit_def)
            print("Hit_accession:", alignment.accession)
            # print("Hit_len:", alignment.len)
            for hsp in alignment.hsps:
                # print("Hsp_bit_score:", hsp.bit_score)
                print("Hsp_score:", hsp.score)
                print("Hsp_evalue:", hsp.expect)
                # print("Hsp_query_from:", hsp.query_from)
                # print("Hsp_query_to:", hsp.query_to)
                # print("Hsp_hit_from:", hsp.hit_from)
                # print("Hsp_hit_to:", hsp.hit_to)
                # print("Hsp_query_frame:", hsp.query_frame)
                # print("Hsp_hit_frame:", hsp.hit_frame)
                # print("Hsp_identity:", hsp.identity)
                # print("Hsp_positive:", hsp.positive)
                print("Hsp_gaps:", hsp.gaps)
                # print("Hsp_align_len:", hsp.align_len)
                # print("Hsp_qseq:", hsp.qseq)
                # print("Hsp_hseq:", hsp.hseq)
                # print("Hsp_midline:", hsp.midline)


def parse_bproc():
    pass


def main():
    run_rpsblast()
    run_rpsbproc()
    # parse_blast()


if __name__ == '__main__':
    main()
