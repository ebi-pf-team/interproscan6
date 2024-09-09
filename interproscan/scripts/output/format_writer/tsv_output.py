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

    with open(output_path, 'w') as tsv_file:
        current_date = datetime.now().strftime('%d-%m-%Y')
        for seq_target, info in seq_matches.items():
            sequence_data = info['sequences']
            matches = info["matches"]
            seq_id = seq_target
            md5 = sequence_data['md5']
            seq_len = sequence_data['length']

            for match_acc, match in matches.items():
                goterms, pathways = [], []
                if match["entry"]:
                    entry_acc = match['entry'].get('accession', "-")
                    sig_desc = match['entry'].get('description') or match['entry'].get('name', "-")
                    entry_desc = match['entry'].get('ipr_description') or match['entry'].get('ipr_name', "-")
                    for go_info in match["entry"]["goXRefs"]:
                        goterms.append(go_info["id"])
                    for pwy_info in match["entry"]["pathwayXRefs"]:
                        pathways.append(pwy_info["id"])
                match_db = match["member_db"]
                xrefs = f"{'|'.join(goterms)}\t{'|'.join(pathways)}"

                for location in match["locations"]:
                    if match["member_db"].upper() in ["SIGNALP", "SIGNALP_EUK"]:
                        sig_acc, status = "Signal Peptide", ""
                        ali_from = match["locations"][0]["start"]
                        ali_to = match["locations"][0]["end"]
                        evalue = match["locations"][0]["pvalue"]
                    elif match_db.upper() in ["CDD", "HAMAP", "PROSITE_PROFILES"]:
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = location["score"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                    elif match_db.upper() in ["COILS", "MOBIDB", "PROSITE_PATTERNS"]:
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = "-"
                        ali_from = location["start"]
                        ali_to = location["end"]
                    elif match_db.upper() == "PIRSF":
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = location["evalue"]
                        ali_from = location["envelopeStart"]
                        ali_to = location["envelopeEnd"]
                    elif match_db.upper() == "PHOBIUS":
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = "-"
                        ali_from = location["start"]
                        ali_to = location["end"]
                        entry_desc = match["description"]
                        if seq_len == ali_from + ali_to - 1:
                            break
                    elif match_db.upper() == "PRINTS":
                        sig_acc = match["accession"]
                        status = "T"
                        evalue = match["evalue"]
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
                        sig_acc, sig_desc, ali_from, ali_to,
                        evalue, status, current_date, entry_acc,
                        entry_desc, xrefs)


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

    with open(output_path, 'w') as tsv_file:
        for seq_target, info in seq_matches.items():
            matches = info["matches"]
            seq_id = seq_target

            for match_acc, match in matches.items():
                member_db = match["member_db"]
                version = match['version']
                try:
                    version_major, version_minor = version.split('.')
                except ValueError:
                    version_major = version
                    version_minor = "0"

                try:  # cdd does not have evalue and score on this level
                    evalue = match["evalue"]
                    score = match["score"]
                except KeyError:
                    evalue = "-"
                    score = "-"

                if 'model-ac' in match:
                    model_ac = match['model-ac']
                elif member_db.upper() in ["SIGNALP", "SIGNALP_EUK"]:  # will probably apply to TMHMM and Phobius when added
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
                    elif member_db.upper() in [
                        "COILS", "MOBIDB", "PHOBIUS",
                        "PROSITE_PATTERNS", "SUPERFAMILY"
                    ]:
                        hmm_start = "-"
                        hmm_end = "-"
                        hmm_length = "-"
                        location_score = "-"
                        env_end, env_start = "-", "-"
                    elif member_db.upper() in ["PRINTS", "SIGNALP", "SIGNALP_EUK"]:
                        hmm_start = location["start"]
                        hmm_end = location["end"]
                        hmm_length = int(hmm_end) - int(hmm_start)
                        location_score = location["pvalue"]
                        env_end, env_start = "-", "-"
                    elif member_db.upper() == "SMART":
                        hmm_start = location["hmmStart"]
                        hmm_end = location["hmmEnd"]
                        hmm_length = location["hmmLength"]
                        location_score = location["score"]
                        env_end = "-"
                        env_start = "-"
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
                    elif member_db.upper() in ["HAMAP", "PROSITE_PROFILES"]:
                        sig_acc = match["accession"]
                        evalue = location["score"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = location["score"]
                    elif member_db.upper() in ["COILS", "MOBIDB", "PHOBIUS", "PROSITE_PATTERNS", "SUPERFAMILY"]:
                        sig_acc = match["accession"]
                        evalue = "-"
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = "-"
                    elif member_db.upper() == "PRINTS":
                        sig_acc = match["accession"]
                        evalue = match["evalue"]
                        ali_from = location["start"]
                        ali_to = location["end"]
                        location_evalue = '-'
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
