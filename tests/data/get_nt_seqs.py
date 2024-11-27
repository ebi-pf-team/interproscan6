#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import shutil
import sys

from pathlib import Path

from tqdm import tqdm

from tests.unit_tests.test_inputs.xrefs.build_entries_json import entries


FASTA = "test_prot.fa"
NT_FASTA = "test_nt.fna"
DL_DIR = Path("temp")
BASE_UNIPROT_URL = "https://rest.uniprot.org/uniparc/search?query=(uniprotkb:<ACC>)"
EMBL_BASE_URL = "https://www.ebi.ac.uk/ena/browser/api/fasta/<DNAID>?download=true"
UNIPROT_ERR_FILE = "failed_uniprot_connections"
DOWNLOAD_ERR_FILE = "failed_fasta_downloads"


def main():
    try:
        DL_DIR.mkdir(parnets=True, exists_ok=True)
    except PermissionError:
        print(f"Do not have permission to make the temp download dir at {DL_DIR}")
        sys.exit(1)

    protein_ids = get_protein_ids(FASTA)
    uniparc_to_gene_ids = get_dna_ids(protein_ids)
    all_dl_paths = download_dna_fasta(uniparc_to_gene_ids)
    concatenate_fasta_files(all_dl_paths)


def get_protein_ids(fasta: str) -> list[str]:
    """Load protein ids from the test_prot.fa file"""
    protein_ids = set()
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):
                protein_ids.add(line.split("|")[1])
    print(f"Loaded {len(protein_ids)} protein ids from {fasta}")
    return list(protein_ids)


def get_dna_ids(protein_ids: list[str]) -> dict[str: set[str]]:
    """Retrieve EMBL DNA uniparc cross references from UniProt"""
    failed_connections = set()
    uniparc_to_gene_ids = {}  # {protein id: set(uniparcCrossReference ids)}
    for prot_id in tqdm(protein_ids, "Checking responses"):
        uniparc_to_gene_ids[prot_id] = set()
        response = requests.get(BASE_UNIPROT_URL.replace("<ACC>", prot_id))
        if response.status_code == 200:
            json_data = response.json()
            for result in json_data["results"]:
                if "uniParcCrossReferences" in result:
                    for result_dict in result["uniParcCrossReferences"]:
                        if result_dict["database"] == "EMBL":
                            uniparc_to_gene_ids[prot_id].add(result_dict["id"])
        else:
            failed_connections.add(prot_id)

    if failed_connections:
        print(
            f"Could not retrieve a response for {len(failed_connections)} protein ids.\n"
            f"These protein ids are logged in {UNIPROT_ERR_FILE}"
        )
        with open(UNIPROT_ERR_FILE, "w") as fh:
            for prot_id in failed_connections:
                fh.write(f"{prot_id}\n")

    return uniparc_to_gene_ids


def download_dna_fasta(uniparc_to_gene_ids: dict[str, set[str]]) -> list[Path]:
    """Download nt FASTA files from EBI/ena"""
    all_dl_paths = []
    failed_connections = []

    for prot_id in tqdm(uniparc_to_gene_ids, desc="Download nt seq fasta for protein ids"):
        for dna_id in uniparc_to_gene_ids[prot_id]:
            url = EMBL_BASE_URL.replace("<DNAID>", dna_id)

            dl_path = DL_DIR / f"{prot_id}_{dna_id}.fna"

            response = requests.get(url)

            if response.status_code == 200:
                all_dl_paths.append(dl_path)
                with open(dl_path, 'wb') as file:
                    file.write(response.content)
            else:
                failed_connections.append(dl_path)

    if all_dl_paths:
        print(f"Successfully downloaded {len(all_dl_paths)} fasta files (available in {DL_DIR})")
    if failed_connections:
        print(
            f"Failed to download {len(failed_connections)} fasta files"
            f"The paths for these files (including the protein id and dna id) are logged in {DOWNLOAD_ERR_FILE}"
        )
        with open(DOWNLOAD_ERR_FILE, "w") as fh:
            for dl_path in failed_connections:
                fh.write(f"{dl_path}\n")

    return all_dl_paths


def concatenate_fasta_files(all_dl_paths: list[Path]) -> None:
    """Concatenate the many download FASTA files into a single file and delete the temp dir"""
    with open(NT_FASTA, "w") as outfile:
        for fasta in all_dl_paths:
            with open(fasta, "r") as infile:
                outfile.write(infile.read())

    print(f"Files have been concatenated into {NT_FASTA}.\nDeleting the temp dir {DL_DIR}")
    shutil.rmtree(DL_DIR)


if __name__ == "__main__":
    main()
