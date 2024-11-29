import json
import os
import shutil
import zipfile
import subprocess
from collections import OrderedDict
from pathlib import Path


def get_current_output(test_output_dir: str,
                       current_output_path: str,
                       input_path: str,
                       applications: str,
                       data_dir: str,
                       disable_precalc: bool) -> dict:
    disable_precalc = "--disable_precalc" if disable_precalc else ""
    command = f"nextflow run main.nf --input {input_path} --applications {applications} " \
              f"--formats json --outdir {test_output_dir} -profile docker \
              --datadir {data_dir}"
    if os.path.exists(str(current_output_path) + ".json"):
        os.remove(str(current_output_path) + ".json")
    subprocess.run(command, shell=True)
    with open(str(current_output_path) + ".json", 'r') as f:
        return json.load(f)


def get_expected_result(expected_result_path: str) -> dict:
    path = Path(expected_result_path)
    if path.suffix == ".zip":
        with zipfile.ZipFile(path, 'r') as zip_ref:
            extracted_dir = path.parent / "extracted_temp"
            zip_ref.extractall(extracted_dir)
            json_files = [f for f in os.listdir(extracted_dir) if f.endswith(".json")]
            json_file = json_files[0]
        expected_result_file = os.path.join(extracted_dir, json_file)
    else:
        expected_result_file = expected_result_path + ".json"
    with open(expected_result_file, 'r') as f:
        expected_result = json.load(f)
    shutil.rmtree(extracted_dir, ignore_errors=True)
    return expected_result


def json2dict(data, ignore_fields):
    match_single_fields = ["modelAccession", "sequenceLength", "evalue", "score", "bias", "graphScan"]
    signature_single_fields = ["accession", "name", "description"]
    entry_fields = ["accession", "name", "description", "type", "goXRefs", "pathwayXRefs"]
    locations_single_fields = ["start", "end", "hmmStart", "hmmEnd", "hmmLength", "hmmBounds", "envelopeStart",
                               "envelopeEnd", "evalue", "score", "bias", "queryAlignment", "targetAlignment",
                               "sequenceFeature", "pvalue", "motifNumber", "level", "cigarAlignment", "representative"]

    match_fields_filtered = [item for item in match_single_fields if item not in ignore_fields]
    signature_fields_filtered = [item for item in signature_single_fields if item not in ignore_fields]
    entry_fields_filtered = [item for item in entry_fields if item not in ignore_fields]
    locations_fields_filtered = [item for item in locations_single_fields if item not in ignore_fields]

    site_single_fields = ["description", "numLocations", "label", "group", "hmmStart", "hmmEnd", "start", "end"]
    site_location_fields = ["residue", "start", "end"]
    location_fragment_fields = ["start", "end", "dcStatus", "dc-status"]
    tree_grafter_fields = ["ancestralNodeID", "graftPoint", "subfamilyAccession", "subfamilyName", "subfamilyDescription", "proteinClass"]

    result = {}
    for result_item in data.get("results", []):
        for match in result_item.get("matches", []):
            signature = match.get("signature", {})
            signature_sorted = OrderedDict((k, signature[k]) for k in signature if k in signature_fields_filtered)
            accession = signature.get("accession")
            library = signature.get("signatureLibraryRelease", {}).get("library", "").lower().replace("-", "").replace(" ", "").replace("_", "")
            if not library:
                continue
            if library not in result:
                result[library] = {}

            for xref in result_item.get("xref", []):
                memberDB = xref.get("id")
                if memberDB not in result[library]:
                    result[library][memberDB] = {}
                if accession not in result[library][memberDB]:
                    result[library][memberDB][accession] = {"locations": {}}

                # Match fields
                for field in match_fields_filtered:
                    if field in match:
                        result[library][memberDB][accession][field] = match[field]

                # Tree grafter fields
                if "treeGrafter" in match and "treeGrafter" not in ignore_fields:
                    result[library][memberDB][accession]["treeGrafter"] = OrderedDict(
                        (k, match["treeGrafter"][k]) for k in match["treeGrafter"] if k in tree_grafter_fields
                    )

                # Signature fields
                result[library][memberDB][accession]["signature"] = signature_sorted
                # Entry fields
                if signature["entry"] and "entry" not in ignore_fields:
                    result[library][memberDB][accession]["signature"]["entry"] = {}
                    for field in entry_fields_filtered:
                        if field in signature["entry"]:
                            result[library][memberDB][accession]["signature"]["entry"][field] = signature["entry"][field]

                # Locations fields
                for loc in match.get("locations", []):
                    location_key = f"{loc.get('start')}-{loc.get('end')}"
                    result[library][memberDB][accession]["locations"][location_key] = OrderedDict((k, loc[k]) for k in loc if k in locations_fields_filtered)
                    if "sites" in loc and "sites" not in ignore_fields:
                        result[library][memberDB][accession]["locations"][location_key]["sites"] = []
                        for site in loc["sites"]:
                            site_sorted = OrderedDict((k, site[k]) for k in site if k in site_single_fields)
                            if "locations" in site:
                                site_sorted["locations"] = []
                                for site_location in site["locations"]:
                                    site_location_sorted = OrderedDict(
                                        (k, site_location[k]) for k in site_location if k in site_location_fields
                                    )
                                    site_sorted["locations"].append(site_location_sorted)
                            result[library][memberDB][accession]["locations"][location_key]["sites"].append(site_sorted)
                        result[library][memberDB][accession]["locations"][location_key]["sites"] = sorted(
                            result[library][memberDB][accession]["locations"][location_key]["sites"],
                            key=lambda x: x.get("description", "") or ""
                        )

                    # location-fragments fields
                    if "location-fragments" in loc and "location-fragments" not in ignore_fields:
                        result[library][memberDB][accession]["locations"][location_key]["location-fragments"] = []
                        for location_fragment in loc["location-fragments"]:
                            location_fragment_sorted = OrderedDict(
                                (k, location_fragment[k]) for k in location_fragment if
                                k in location_fragment_fields
                            )
                            result[library][memberDB][accession]["locations"][location_key][
                                "location-fragments"].append(location_fragment_sorted)

                sorted_locations = OrderedDict(
                    sorted(
                        result[library][memberDB][accession]["locations"].items(),
                        key=lambda x: int(x[0].split('-')[0])
                    )
                )
                result[library][memberDB][accession]["locations"] = sorted_locations

    return OrderedDict(sorted(result.items()))


def compare(expected: dict,
            current: dict,
            output_file):
        all_libraries = set(expected.keys()).union(set(current.keys()))
        with open(output_file, "w") as file:
            for library in all_libraries:
                if library not in current:
                    continue
                all_xref_ids = set(expected[library].keys()).union(set(current[library].keys()))
                for xref_id in all_xref_ids:
                    if xref_id not in expected[library]:
                        file.write(f"EXTRA SEQUENCE: {library}: xref.id '{xref_id}'\n")
                        continue
                    if xref_id not in current[library]:
                        file.write(f"MISSING SEQUENCE: {library}: xref.id '{xref_id}'\n")
                        continue
                    all_accessions = set(expected[library][xref_id].keys()).union(set(current[library][xref_id].keys()))
                    for accession in all_accessions:
                        if accession not in expected[library][xref_id]:
                            file.write(f"EXTRA ACCESSION: {library} -> {xref_id}: accession '{accession}'\n")
                            continue
                        if accession not in current[library][xref_id]:
                            file.write(f"MISSING ACCESSION: {library} -> {xref_id}: accession '{accession}'\n")
                            continue
                        for field in expected[library][xref_id][accession]:
                            value1 = expected[library][xref_id][accession][field]
                            value2 = current[library][xref_id][accession].get(field)
                            if value1 != value2:
                                file.write(f"MISMATCH ON {field.upper()}: {library} -> {xref_id} -> {accession}:\n\t Expected: {value1} \n\t Current : {value2}\n")


def test_json_output(test_output_dir, input_path, expected_result_path, output_path, applications, data_dir, disable_precalc):
    expected_result = get_expected_result(expected_result_path)
    current_result = get_current_output(test_output_dir, output_path, input_path, applications, data_dir, disable_precalc)

    ignore_fields = ["description", "name", "representative", "evalue"]
    expected = json2dict(expected_result, ignore_fields)
    current = json2dict(current_result, ignore_fields)

    output_file = "tests/integration_tests/mismatches.txt"
    compare(expected, current, output_file)

    assert os.stat(output_file).st_size == 0
