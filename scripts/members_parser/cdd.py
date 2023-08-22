from Bio.Blast.Applications import NcbirpsblastCommandline
import subprocess
import re


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


def parse_rpsbproc(rpsbproc_file):
    # DATA_BLOCK_START_MARKER = "DATA"
    # DATA_BLOCK_END_MARKER = "ENDDATA"
    SESSION_BLOCK_START_MARKER = "SESSION"
    SESSION_BLOCK_END_MARKER = "ENDSESSION"
    QUERY_BLOCK_START_MARKER = "QUERY"
    # QUERY_BLOCK_END_MARKER = "ENDQUERY"

    DOMAINS_BLOCK_START_MARKER = "DOMAINS"
    DOMAINS_BLOCK_END_MARKER = "ENDDOMAINS"
    SITES_BLOCK_START_MARKER = "SITES"
    SITES_BLOCK_END_MARKER = "ENDSITES"

    # MOTIFS_BLOCK_START_MARKER = "MOTIFS"
    # MOTIFS_BLOCK_END_MARKER = "ENDMOTIFS"

    QUERY_LINE_PATTERN = re.compile("^QUERY\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(.*)$")
    DOMAIN_LINE_PATTERN = re.compile("^(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")
    SITE_LINE_PATTERN = re.compile("^(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\d+)\\s+(\\d+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)")

    with open(rpsbproc_file, 'r') as f:
        lines = f.readlines()

    all_matches = {}
    matches = {}
    sites = {}
    protein_identifier = ""
    # pssmid2modelId = {}
    definitionLine = ""

    domain_block_state = False
    site_block_state = False
    for line in lines:
        line = line.strip()

        if line.startswith("#"):
            pass
        else:
            if line.startswith(SESSION_BLOCK_START_MARKER):
                protein_identifier = line.split()[1].strip()

            elif line.startswith(SESSION_BLOCK_END_MARKER):
                print(protein_identifier)
                all_matches[protein_identifier] = matches
                all_matches[protein_identifier] = sites

            elif line.startswith(QUERY_BLOCK_START_MARKER):
                match = QUERY_LINE_PATTERN.match(line)
                if match:
                    matches["queryId"] = match.group(1)
                    matches["sequenceType"] = match.group(2)
                    matches["sequenceLength"] = int(match.group(3))
                    matches["definitionLine"] = match.group(4)
                    matches["sequenceIdentifier"] = definitionLine.strip()

            elif line.startswith(DOMAINS_BLOCK_START_MARKER) or domain_block_state:
                domain_block_state = True
                match = DOMAIN_LINE_PATTERN.match(line)
                if match:
                    matches["sessionNumber"] = int(match.group(1))
                    matches["queryId"] = match.group(2)
                    matches["queryType"] = match.group(3)
                    # match["hitType"] = queryType  # !!
                    matches["pssmID"] = match.group(4)
                    matches["locationStart"] = int(match.group(5))
                    matches["locationEnd"] = int(match.group(6))
                    matches["eValue"] = float(match.group(7))
                    matches["score"] = float(match.group(8))
                    matches["model"] = match.group(9)
                    matches["shortName"] = match.group(10)
                    matches["incomplete"] = match.group(11)
                    matches["superfamilyPSSMId"] = match.group(12)
                    # match["pssmid2modelId"][pssmID] = model
            elif line.startswith(DOMAINS_BLOCK_END_MARKER):
                domain_block_state = False

            elif line.startswith(SITES_BLOCK_START_MARKER) or site_block_state:
                site_block_state = True
                siteInfo = line.split('\t')
                if len(siteInfo) > 5:
                    sites["sessionNumber"] = int(siteInfo[0])
                    sites["queryId"] = siteInfo[1]
                    sites["annotQueryType"] = siteInfo[2]
                    # sites["annotationType"] = annotQueryType  # !!
                    sites["title"] = siteInfo[3]
                    sites["residues"] = siteInfo[4]
                    sites["completeSize"] = int(siteInfo[5])
                    sites["mappedSize"] = int(siteInfo[6])
                    sites["sourceDomain"] = siteInfo[7]
                    # sites["model"] = pssmid2modelId.get(sourceDomain)

            elif line.startswith(SITES_BLOCK_END_MARKER):
                    site_block_state = False


def main():
    run_rpsblast(all_appl, output_blast_path)
    run_rpsbproc(output_blast_path, output_rpsbproc_path)
    parse_rpsbproc(output_rpsbproc_path)


if __name__ == '__main__':
    main()
