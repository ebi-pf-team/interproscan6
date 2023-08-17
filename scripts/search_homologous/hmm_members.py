import pyhmmer
from pyhmmer.easel import SequenceFile, Alphabet
from pyhmmer.plan7 import HMMFile
from Bio import SeqIO
import hashlib
import argparse

# hmmsearch search profile against sequence database
# hmmscan search sequence against profile database
# hmmpress prepare profile database for hmmscan

# [AntiFam-7.0,CDD-3.20,Coils-2.2.1,FunFam-4.3.0,Gene3D-4.3.0,Hamap-2023_01,MobiDBLite-2.0,NCBIfam-12.0,PANTHER-17.0,Pfam-35.0,PIRSF-3.10,PIRSR-2021_05,PRINTS-42.0,ProSitePatterns-2022_05,ProSiteProfiles-2022_05,SFLD-4,SMART-9.0,SUPERFAMILY-1.75]
# hmmsearch -Z 61295632 --cut_ga --cpu 1 -o output_hmmer.json ./data/pfam/35.0/pfam_a.search_homologous ./files_test/test_all_appl.fasta

# WHEN TESTS WITH THIS MEMBER FINISH THIS PARAMETERS WILL GO TO nextflow.config
project_dir = "/Users/lcf/interproscan6"
all_appl = f"{project_dir}/files_test/test_all_appl.fasta"

############################
# Member databases HMMs
############################
# Families, domains, sites & repeats
hmm_hamap_path = f"{project_dir}/data/hamap/2023_01/hamap.hmm.lib"
hmm_ncbifam_path = f"{project_dir}/data/ncbifam/11.0/ncbifam.hmm"
hmm_panther_path = f"{project_dir}/data/panther/17.0/famhmm/binHmm"
hmm_pfam_path = f"{project_dir}/data/pfam/35.0/pfam_a.hmm"
hmm_pirsf_path = f"{project_dir}/data/pirsf/3.10/sf_hmm_all"
hmm_sfld_path = f"{project_dir}/data/sfld/4/sfld.hmm"
# Structural domains
hmm_gene3d_path = f"{project_dir}/data/gene3d/4.3.0/gene3d_main.hmm"
hmm_superfamily_path = f"{project_dir}/data/superfamily/1.75/hmmlib_1.75"

# Other category
hmm_antifam_path = f"{project_dir}/data/antifam/7.0/AntiFam.hmm"
hmm_pirsr_path = f"{project_dir}/data/pirsr/2021_05/sr_hmm_all"

# Families, domains, sites & repeats
hmm_prositeprofiles_path = ""
hmm_prositepatterns_path = ""
hmm_smart_path = f"{project_dir}/data/smart/9.0/smart.HMMs"
hmm_prints_path = ""

# Other sequence features
#     phobius, signalp, coils, mobidblite, tmhmm
# Other category
hmm_signalp_euk_path = ""
hmm_signalp_gram_positive_path = ""
hmm_signalp_gram_negatove_path = ""
hmm_funfam_file_10 = f"{project_dir}/data/funfam/4.3.0/models/4/10/75/10.hmm"
hmm_funfam_file_40 = f"{project_dir}/data/funfam/4.3.0/models/3/40/50/40.hmm"
hmm_funfam_file_1170 = f"{project_dir}/data/funfam/4.3.0/models/3/40/50/1170.hmm"

hmm_paths = {
    "pfam": hmm_pfam_path,
    "antifam": hmm_antifam_path,
    "ncbifam": hmm_ncbifam_path,
    "hamap": hmm_hamap_path,
    "pirsr": hmm_pirsr_path,
    "pirsf": hmm_pirsf_path,
}


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


def search_matches(fasta_file, member_list: set):
    alphabet = Alphabet.amino()
    # background = pyhmmer.plan7.Background(alphabet)
    # pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

    with SequenceFile(fasta_file, digital=True, alphabet=alphabet) as sf:
        sequencesD = sf.read_block()

    homologs_result = {}
    for appl in member_list:
        appl_result = {}
        with HMMFile(hmm_paths[appl]) as f:
            hmm_list = list(f)
        for sequence in pyhmmer.hmmsearch(hmm_list, sequencesD):
            sequence_results = {}
            hits_result = {}
            for hit in sequence:
                hit_name = hit.name.decode()
                if hit_name in hits_result:
                    hits_result[hit_name].append(format_hit(hit))
                else:
                    hits_result[hit_name] = [format_hit(hit)]

                domains_result = {}
                for domain in hit.domains:
                    if domains_result:
                        domains_result.append(format_domain(domain))
                    else:
                        domains_result = [format_domain(domain)]

                hits_result["domains"] = domains_result

            sequence_results["hits"] = hits_result

            sequence_id = sequence.query_name.decode()
            appl_result[sequence_id] = sequence_results
        homologs_result[appl] = appl_result
    return homologs_result


def format_seq_result(hits):
    formatted_hits = {
        "query_name": hits.query_name.decode(),
        "query_accession": hits.query_accession,
        # "reported": hits.reported,
        # "included": hits.included,
        "block_length": hits.block_length,
        # "domE": hits.domE,
        # "domT": hits.domT,
        # "domZ": hits.domZ,
        # "incdomE": hits.incdomE,
        # "incdomT": hits.incdomT,
        "bit_cutoffs": hits.bit_cutoffs,
        "searched_models": hits.searched_models,
        "searched_nodes": hits.searched_nodes,
        "searched_residues": hits.searched_residues,
        "searched_sequences": hits.searched_sequences,
        "long_targets": hits.long_targets,
        "strand": hits.strand,
        # "Z": hits.Z,
        # "T": hits.T,
        # "E": hits.E,
        # "incE": hits.incE,
        # "incT": hits.incT
    }
    return formatted_hits


def format_hit(hit):
    formatted_hit = {
        "Hit name": hit.name.decode(),
        "accession": hit.accession,
        "description": hit.description,
        "E-value": hit.evalue,
        "Score": hit.score,
        "Bias": hit.bias,
        "Reported": hit.reported,
        "dropped": hit.dropped,
        "duplicate": hit.duplicate,
        "included": hit.included,
        "new": hit.new,
        "pvalue": hit.pvalue,
        "pre_score": hit.pre_score,
        "sum_score": hit.sum_score,
        "number_of_domains_found": len(hit.domains),
        "best_domain_bias": hit.best_domain.bias,
        "best_domain_score": hit.best_domain.score,
        "best_domain_pvalue": hit.best_domain.pvalue,
        "best_domain_included": hit.best_domain.included,
        "best_domain_reported": hit.best_domain.reported,
        "best_domain_c_evalue": hit.best_domain.c_evalue,
        "best_domain_correction": hit.best_domain.correction,
        "best_domain_env_from": hit.best_domain.env_from,
        "best_domain_env_to": hit.best_domain.env_to,
        "best_domain_envelope_score": hit.best_domain.envelope_score,
        "best_domain_i_evalue": hit.best_domain.i_evalue,
    }
    return formatted_hit


def format_domain(domain):
    domain_formated = {
        "bias": domain.bias,
        "score": domain.score,
        "pvalue": domain.pvalue,
        "reported": domain.reported,
        "correction": domain.correction,
        "c_evalue": domain.c_evalue,
        "env_from": domain.env_from,
        "env_to": domain.env_to,
        "envelope_score": domain.envelope_score,
        "included": domain.included,
        "i_evalue": domain.i_evalue,
        # "domain_coverage": (hit.best_domain.env_to - hit.best_domain.env_from + 1) / model_length

    }
    return domain_formated


def main():
    members_with_hmm = ["pfam", "antifam", "ncbifam", "pirsf", "pirsr", "hamap"]

    parser = argparse.ArgumentParser(
        description="Request to calc matches for members with hmm"
    )
    parser.add_argument("-fasta", "--fastafile", type=str, help="fasta file with sequences")
    parser.add_argument("-appl", "--applications", nargs="*", help="list of analysis")
    args = parser.parse_args()

    applications = args.applications
    applications_hmm = set(applications) & set(members_with_hmm)

    sequences = get_sequences(all_appl)
    homologs_found = search_matches(all_appl, applications_hmm)
    print(sequences)
    print(homologs_found)


if __name__ == "__main__":
    main()
