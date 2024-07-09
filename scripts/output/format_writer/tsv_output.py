import os
from datetime import datetime


def tsv_output(seq_matches: dict, output_path: str):
    def write_to_tsv(
            seq_id, md5, seq_len, member_db, sig_acc,
            sig_desc, ali_from, ali_to, evalue, status,
            current_date, interpro_acc, interpro_name, xrefs):
        tsv_file.write((
            f"{seq_id}\t{md5}\t{seq_len}\t{member_db}\t{sig_acc}\t"
            f"{sig_desc}\t{ali_from}\t{ali_to}\t{evalue}\t{status}\t"
            f"{current_date}\t{interpro_acc}\t{interpro_name}\t{xrefs}\n"
        ))

    tsv_output = os.path.join(output_path + '.tsv')

    with open(tsv_output, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        for seq_target, info in seq_matches.items():
            sequence_data = info['sequences']
            matches = info["matches"]
            seq_id = seq_target
            md5 = sequence_data[2]
            seq_len = sequence_data[3]

            for match_acc, match in matches.items():
                entry_acc, entry_name, entry_desc = "-", "-", "-"
                goterms, pathways = [], []
                if match["entry"]:
                    entry_acc = match["entry"]["accession"]
                    entry_name = match["entry"]["short_name"]
                    entry_desc = match["entry"]["name"]
                    for go_info in match["entry"]["goXRefs"]:
                        goterms.append(go_info["id"])
                    for pwy_info in match["entry"]["pathwayXRefs"]:
                        pathways.append(pwy_info["id"])
                match_db = match["member_db"]
                xrefs = f"{'|'.join(goterms)}\t{'|'.join(pathways)}"

                for location in match["locations"]:
                    if match["member_db"].upper() == "SIGNALP":
                        sig_acc, status = "Signal Peptide", ""
                        ali_from = match["locations"][0]["start"]
                        ali_to = match["locations"][0]["end"]
                        evalue = match["locations"][0]["pvalue"]

                    if match["member_db"].upper() == "DEEPTMHMM":
                        sig_acc, status = location["location_tag"], ""
                        ali_from = location["start"]
                        ali_to = location["end"]
                        evalue = "-"
                        status = "T"

                    elif match_db.upper() in ["CDD", "HAMAP", "PROSITE_PROFILES"]:
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = location["score"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                    elif match_db.upper() == "PROSITE_PATTERNS":
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = "-"
                        ali_from = location["start"]
                        ali_to = location["end"]
                    else:
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]

                    write_to_tsv(
                        seq_id, md5, seq_len, match_db,
                        sig_acc, entry_desc, ali_from, ali_to,
                        evalue, status, current_date, entry_acc,
                        entry_name, xrefs)


def tsv_pro_output(seq_matches: dict, output_path: str):
    def write_to_tsv(
            member_db, version_major, version_minor, seq_id,
            sig_acc, model_ac, ali_from, ali_to, fragment,
            score, evalue, raw_hmm_bound, hmm_start, hmm_end,
            hmm_length, env_start, env_end, location_score,
            location_evalue, cigar_alignment):
        tsv_file.write((
            f"{member_db}\t{version_major}\t{version_minor}\t{seq_id}\t"
            f"{sig_acc}\t{model_ac}\t{ali_from}\t{ali_to}\t{fragment}\t"
            f"{score}\t{evalue}\t{raw_hmm_bound}\t{hmm_start}\t{hmm_end}\t"
            f"{hmm_length}\t{env_start}\t{env_end}\t{location_score}\t"
            f"{location_evalue}\t{cigar_alignment}\n"
        ))

    tsv_output = os.path.join(output_path + '.tsv-pro')

    with open(tsv_output, 'w') as tsv_file:
        for seq_target, info in seq_matches.items():
            matches = info["matches"]
            seq_id = seq_target

            for match_acc, match in matches.items():
                member_db = match["member_db"]
                try:
                    version_major, version_minor = match['version'].split('.')
                except ValueError:
                    version_major = match['version']
                    version_minor = "0"
                try:  # cdd does not have evalue and score on this level
                    evalue = match["evalue"]
                    score = match["score"]
                except KeyError:
                    evalue = "-"
                    score = "-"

                if 'model-ac' in match:
                    model_ac = match['model-ac']
                elif member_db.upper() in ["SIGNALP", "DEEPTMHMM"]:  # will probably apply to TMHMM and Phobius when added
                    model_ac = "-"
                else:
                    model_ac = match['accession']

                for location in match["locations"]:
                    if member_db.upper() in ["CDD", "HAMAP", "PROSITE_PROFILES"]:
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["score"]
                        env_end, env_start = "-", "-"
                    elif member_db.upper() == "PROSITE_PATTERNS":
                        hmm_start = "-"
                        hmm_end = "-"
                        hmm_length = "-"
                        location_score = "-"
                        env_end, env_start = "-", "-"
                    elif member_db.upper() == "SIGNALP":
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["pvalue"]
                        env_end, env_start = "-", "-"
                    elif member_db.upper() == "DEEPTMHMM":
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["location_tag"]
                        env_end, env_start = "-", "-"
                    else:
                        hmm_start = location["hmmStart"]
                        hmm_end = location["hmmEnd"]
                        hmm_length = location["hmmLength"]
                        location_score = location["score"]
                        env_end = location["envelopeEnd"]
                        env_start = location["envelopeStart"]
                    try:
                        fragment = location["fragment"]
                    except KeyError:
                        fragment = f"{location['start']}-{location['end']}-S"
                    try:
                        raw_hmm_bound = location["rawHmmBounds"]
                    except KeyError:
                        raw_hmm_bound = ""  # lookup match does not have rawHmmBounds

                    if match_acc == "signal_peptide":
                        sig_acc, status = "Signal Peptide", ""
                        ali_from = match["locations"][0]["start"]
                        ali_to = match["locations"][0]["end"]
                        location_evalue = match["locations"][0]["pvalue"]
                    elif member_db.upper() == "DEEPTMHMM":
                        sig_acc, status = location["location_tag"], ""
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = "-"
                    elif member_db.upper() in ["HAMAP", "PROSITE_PROFILES"]:
                        sig_acc = match["accession"]
                        evalue = location["score"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = location["score"]
                    elif member_db.upper() == "PROSITE_PATTERNS":
                        sig_acc = match["accession"]
                        evalue = "-"
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = "-"
                    else:
                        sig_acc = match["accession"]
                        evalue = location["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = location["evalue"]
                    cigar_alignment = ""
                    try:
                        cigar_alignment = location["cigar_alignment"]
                    except KeyError:
                        pass  # some members may not have cigar alignment

                    write_to_tsv(
                        member_db, version_major, version_minor, seq_id,
                        sig_acc, model_ac, ali_from, ali_to, fragment,
                        score, evalue, raw_hmm_bound, hmm_start, hmm_end,
                        hmm_length, env_start, env_end, location_score,
                        location_evalue, cigar_alignment)

