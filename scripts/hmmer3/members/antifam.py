import pyhmmer
from pyhmmer.easel import SequenceFile
from pyhmmer.plan7 import HMMFile

hmm_pfam_file = "/Users/lcf/interproscan6/data/pfam/35.0/pfam_a.hmm" #temp! just to test with the hit in this input
hmm_antifam_file = "/Users/lcf/interproscan6/data/antifam/7.0/AntiFam.hmm"
fasta_file = "/Users/lcf/interproscan6/files_test/test_all_appl.fasta"

# <property name="fullPathToHmmScanBinary" value="${binary.hmmer3.hmmscan.path}"/>
# <property name="fullPathToHmmsearchBinary" value="${binary.hmmer33.hmmsearch.path}"/>
# <property name="binarySwitches" value="${hmmer3.hmmsearch.switches.antifam} ${hmmer3.hmmsearch.cpu.switch.antifam}"/>


def search_homologs(fasta_file, hmm_antifam_file):
    alphabet = pyhmmer.easel.Alphabet.amino()
    background = pyhmmer.plan7.Background(alphabet)
    pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

    with SequenceFile(fasta_file, digital=True) as seq_file:
        sequences = seq_file.read_block()
    # for seq in sequences:
    #     print(seq.name)

    with HMMFile(hmm_antifam_file) as hmm_file:
        hmm = hmm_file.read()
        # for hits in pyhmmer.hmmsearch(hmm_file, sequences):
        #     name = hits.query_name.decode()
        #     print(f"HMM {name} found {len(hits)} hits")

    for hits in pipeline.search_hmm(hmm, sequences):
        print(format_hit(hits))

    # hits = pipeline.search_hmm(hmm, seq_file)

    # hits = []
    # hits.append(pipeline.search_hmm(hmm, block))

    # import collections
    # Result = collections.namedtuple("Result", ["query", "cog", "bitscore"])
    #
    # results = []
    # for hits in pyhmmer.hmmsearch(hmm, sequences):
    #     cog = hits.query_name.decode()
    #     for hit in hits:
    #         results.append(Result(hit.name.decode(), cog, hit.score))


def format_hit(hit):
    formatted_hit = {
        "accession": hit.accession,
        "description": hit.description,
        "Hit name": hit.name,
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
        "best_domain_i_evalue": hit.best_domain.i_evalue
    }

        # "start": hit.query_start,
        # "end": hit.query_end,
        # "coverage": (hit.query_end - hit.query_start + 1) / hit.model_length

        # "query_name": hit.query_name,
        # "query_accession": hit.query_accession,
        # "reported": hit.reported,
        # "included": hit.included,
        # "block_length": hit.block_length,
        # "domE": hit.domE,
        # "domT": hit.domT,
        # "domZ": hit.domZ,
        # "incdomE": hit.incdomE,
        # "incdomT": hit.incdomT,
        # "bit_cutoffs": hit.bit_cutoffs,
        # "searched_models": hit.searched_models,
        # "searched_nodes": hit.searched_nodes,
        # "searched_residues": hit.searched_residues,
        # "searched_sequences": hit.searched_sequences,
        # "long_targets": hit.long_targets,
        # "strand": hit.strand,
        # "Z": hit.Z,
        # "T": hit.T,
        # "E": hit.E,
        # "incE": hit.incE,
        # "incT": hit.incT

    return formatted_hit


if __name__ == '__main__':
    search_homologs(fasta_file, hmm_pfam_file)
