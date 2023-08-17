from Bio.Blast.Applications import NcbirpsblastCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Blast import NCBIXML
import subprocess


# WHEN TESTS WITH THIS MEMBER FINISH THIS PARAMETERS WILL GO TO nextflow.config
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
output_rpsbproc_path = f"{project_dir}/results/cdd_blast_output/output_rpsbproc.txt"


def run_rpsblast(sequence: str, output_path: str):
    cmd = NcbirpsblastCommandline(
        cmd=binary_rpsblast_path,
        query=sequence,
        db=cdd_library_path,
        evalue=0.01,
        seg="no",
        outfmt=11,
        out=output_path
    )
    cmd()


def run_rpsbproc(output_blast: str, output_path: str):
    command = [
        binary_rpsbproc_path,
        "-i", output_blast,
        "-d", cdd_data_path, "--data-mode", "std"]

    with open(output_path, "w") as out:
        try:
            subprocess.run(command, check=True, stdout=out)
        except subprocess.CalledProcessError as e:
            print(e)


def parse_results(output_rpsbproc: str):
    result_handle = open(output_rpsbproc)
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

    # private static final String DATA_BLOCK_START_MARKER = "DATA";
    # private static final String DATA_BLOCK_END_MARKER = "ENDDATA";
    # private static final String SESSION_BLOCK_START_MARKER = "SESSION";
    # private static final String SESSION_BLOCK_END_MARKER = "ENDSESSION";
    # private static final String QUERY_BLOCK_START_MARKER = "QUERY";
    # private static final String QUERY_BLOCK_END_MARKER = "ENDQUERY";
    #
    # private static final String DOMAINS_BLOCK_START_MARKER = "DOMAINS";
    # private static final String DOMAINS_BLOCK_END_MARKER = "ENDDOMAINS";
    # private static final String SITES_BLOCK_START_MARKER = "SITES";
    # private static final String SITES_BLOCK_END_MARKER = "ENDSITES";
    #
    # private static final String MOTIFS_BLOCK_START_MARKER = "MOTIFS";
    # private static final String MOTIFS_BLOCK_END_MARKER = "ENDMOTIFS";
    #
    #     /**
    #      * #QUERY	<query-id>	<seq-type>	<seq-length>	<definition-line>
    #      * #DOMAINS
    #      * #<session-ordinal>	<query-id[readingframe]>	<hit-type>	<PSSM-ID>	<from>	<to>	<E-Value>	<bitscore>	<accession>	<short-name>	<incomplete>	<superfamily PSSM-ID>
    #      * QUERY	Query_1	Peptide	590	sp|Q96N58|ZN578_HUMAN Zinc finger protein 578 OS=Homo sapiens GN=ZNF578 PE=2 SV=2
    #      * DOMAINS
    #      * 1	Query_1	Specific	143639	24	60	3.46102e-15	69.5006	cd07765	KRAB_A-box	-	271597
    #      * ENDDOMAINS
    #      *
    #      *SITES
    #      #<session-ordinal>	<query-id[readingframe]>	<annot-type>	<title>	<residue(coordinates)>	<complete-size>	<mapped-size>	<source-domain>
    #      1	Query_1	Specific	ATP binding site	P272,P273,G274,T275,G276,K277,T278,L279,D330,N377	10	10	99707
    #      1	Query_1	Specific	Walker A motif	G271,P272,P273,G274,T275,G276,K277,T278	8	8	99707
    #      1	Query_1	Specific	arginine finger	R885	1	1	99707
    #      ENDSITES
    #      */
    #
    # private static final Pattern QUERY_LINE_PATTERN
    #         = Pattern.compile("^QUERY\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(.*)$");
    # private static final Pattern DOMAIN_LINE_PATTERN
    #         = Pattern.compile("^(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");
    #
    # private static final Pattern SITE_LINE_PATTERN
    #         =  Pattern.compile("^(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)");


def main():
    run_rpsblast(all_appl, output_blast_path)
    run_rpsbproc(output_blast_path, output_rpsbproc_path)
    # parse_results(output_rpsbproc_path)


if __name__ == '__main__':
    main()
