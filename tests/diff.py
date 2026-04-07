#!/usr/bin/env python

import argparse
import gzip
import json
import shutil
import xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Any

Hit = tuple[
    str,  # Seq ID
    str,  # Souce
    str,  # Signature database[/version]
    str,  # Signature accession
    str,  # InterPro accession
    str,  # comma-separated, ordered, GO terms
    str,  # comma-separated locations or fragments
]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--representative-locations-only", action="store_true")
    parser.add_argument("expected", type=Path)
    parser.add_argument("observed", type=Path)

    args = parser.parse_args()

    with TemporaryDirectory() as tmpdir:
        exp_file = decompress_if_gzipped(args.expected, tmpdir)
        obs_file = decompress_if_gzipped(args.observed, tmpdir)

        extended = supports_extended(exp_file) and supports_extended(obs_file)
        repr_locs_only = args.representative_locations_only

        exp_hits = parse(exp_file, extended=extended, repr_locs_only=repr_locs_only)
        obs_hits = parse(obs_file, extended=extended, repr_locs_only=repr_locs_only)

    print(len(obs_hits & exp_hits), "/", len(exp_hits))
    assert exp_hits == obs_hits


def supports_extended(file: Path) -> bool:
    return file.name.endswith((".json", ".jsonl", ".xml"))


def decompress_if_gzipped(file: Path, tmpdir: str) -> Path:
    if file.suffix != ".gz":
        return file

    tmpfile = Path(tmpdir) / file.with_suffix("").name
    with gzip.open(file, "rb") as src, tmpfile.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)

    return tmpfile


def parse(file: Path, extended: bool, repr_locs_only: bool) -> set[Hit]:
    fmt = file.suffix[1:]
    if fmt == "gff3":
        return parse_gff3(file, repr_locs_only=repr_locs_only)
    elif fmt == "json":
        return parse_json(file, extended=extended, repr_locs_only=repr_locs_only)
    elif fmt == "jsonl":
        return parse_jsonl(file, extended=extended, repr_locs_only=repr_locs_only)
    elif fmt == "tsv":
        return format_tsv(file)
    elif fmt == "xml":
        return format_xml(file, extended=extended)
    else:
        raise NotImplementedError(f"Format not supported: {fmt}")


def parse_gff3(file: Path, repr_locs_only: bool) -> set[Hit]:
    hits = set()

    with file.open("rt") as fh:
        for line in map(str.rstrip, fh):
            if line.startswith("##"):
                if line == "##FASTA":
                    break

                continue

            values = line.split("\t")
            assert len(values) == 9

            seq_id = values[0]
            source = values[1]
            start = int(values[3])
            end = int(values[4])
            attrs = values[8]
            if source == "esl-translate":
                # TODO: support GFF3 parsing for nucleic sequences
                continue

            d_attrs = defaultdict(list)
            for attr in attrs.split(";"):
                key, value = attr.split("=")
                d_attrs[key].append(value)

            if repr_locs_only and d_attrs["representative"][0] == "false":
                continue

            siglib = sig_acc = ""
            if source == "InterPro-N":
                for dbxref in d_attrs["Dbxref"]:
                    dbkey, dbval = dbxref.split(":", maxsplit=1)
                    if dbkey != "InterPro":
                        siglib = dbkey
                        sig_acc = dbval
                        break
            else:
                siglib = source
                sig_acc = d_attrs.get("Alias", d_attrs["Name"])[0]

            interpro_acc = ""
            for dbxref in d_attrs["Dbxref"]:
                dbkey, dbval = dbxref.split(":", maxsplit=1)
                if dbkey == "InterPro":
                    interpro_acc = dbval
                    break

            if d_attrs.get("Ontology_term"):
                go_terms = d_attrs["Ontology_term"][0].split(",")
            else:
                go_terms = []

            hits.add(
                (
                    seq_id,
                    source,
                    siglib,
                    sig_acc,
                    interpro_acc,
                    ",".join(sorted(go_terms)),
                    f"{start}-{end}",
                )
            )

    return hits


def format_tsv(file: str) -> set[Hit]:
    hits = set()
    with file.open("rt") as fh:
        for line in map(str.rstrip, fh):
            values = line.split("\t")
            assert len(values) == 15
            seq_id = values[0]
            sig_lib = values[3]
            sig_acc = values[4]
            start = int(values[6])
            end = values[7]
            source = values[9]
            interpro_acc = values[11] if values[11] != "-" else ""
            go_terms = set()
            if values[13] != "-":
                for go_term in values[13].split("|"):
                    go_id, _ = go_term.split("(")  # GO:0003700(InterPro)
                    go_terms.add(go_id)

            hits.add((
                seq_id, 
                source, 
                sig_lib,
                sig_acc, 
                interpro_acc,
                ",".join(sorted(go_terms)),
                f"{start}-{end}"
            ))

    return hits


def parse_json(file: Path, extended: bool, repr_locs_only: bool) -> set[Hit]:
    with file.open("rt") as fh:
        data = json.load(fh)

    if "results" in data:
        return parse_external_json(
            data["results"], extended=extended, repr_locs_only=repr_locs_only
        )
    else:
        return parse_internal_json(
            data, extended=extended, repr_locs_only=repr_locs_only
        )


def parse_external_json(
    data: list[dict[str, Any]], extended: bool, repr_locs_only: bool
) -> set[tuple[str, str, str, str]]:
    hits = set()
    for seq in data:
        for xref in seq["xref"]:
            seq_id = xref["id"]
            for match in seq["matches"]:
                libojb = match["signature"]["signatureLibraryRelease"]
                siglib = libojb["library"] or ""
                libver = libojb["version"] or ""
                interpro_acc = ""
                go_terms = set()

                if match["signature"]["entry"]:
                    interpro_acc = match["signature"]["entry"]["accession"]
                    for go_term in match["signature"]["entry"]["goXRefs"]:
                        go_terms.add(go_term["id"])

                if siglib == "PANTHER":
                    for go_term in match.get("goXRefs", []):
                        go_terms.add(go_term["id"])

                for loc in match["locations"]:
                    fragments = []
                    if extended:
                        for f in sorted(
                            loc["location-fragments"],
                            key=lambda x: (x["start"], x["end"]),
                        ):
                            values = []
                            for key in ["start", "end", "dc-status"]:
                                value = f[key]
                                values.append(str(value))

                            fragments.append("-".join(values))
                    else:
                        start = loc["start"]
                        end = loc["end"]
                        fragments.append(f"{start}-{end}")

                    hits.add(
                        (
                            seq_id,
                            match["source"],
                            f"{siglib}/{libver}" if extended else siglib,
                            match["signature"]["accession"],
                            interpro_acc,
                            ",".join(sorted(go_terms)),
                            ",".join(fragments),
                        )
                    )

    return hits


def parse_internal_json(
    data: dict[str, dict[str, dict[str, Any]]], extended: bool, repr_locs_only: bool
) -> set[Hit]:
    hits = set()
    for seq_id, seq_matches in data.items():
        if not seq_matches:
            # No matches found in this sequence
            hits.add((seq_id, "", "", "", "", "", ""))
            continue

        for match in seq_matches.values():
            libojb = match["signature"]["signatureLibraryRelease"]
            siglib = libojb["library"] or ""
            libver = libojb["version"] or ""
            interpro_acc = ""
            go_terms = set()

            if match["signature"]["entry"]:
                interpro_acc = match["signature"]["entry"]["accession"]
                for go_term in match["signature"]["entry"]["goXRefs"]:
                    go_terms.add(go_term["id"])

            if match["treegrafter"]:
                for go_term in match["treegrafter"]["goXRefs"]:
                    go_terms.add(go_term["id"])

            for loc in match["locations"]:
                if repr_locs_only and not loc["representative"]:
                    continue

                fragments = []
                if extended:
                    for f in sorted(
                        loc["fragments"], key=lambda x: (x["start"], x["end"])
                    ):
                        values = []
                        for key in ["start", "end", "dcStatus"]:
                            value = f[key]
                            values.append(str(value))

                        fragments.append("-".join(values))
                else:
                    start = loc["start"]
                    end = loc["end"]
                    fragments.append(f"{start}-{end}")

                hits.add(
                    (
                        seq_id,
                        match["source"],
                        match["signature"]["accession"],
                        f"{siglib}/{libver}" if extended else siglib,
                        interpro_acc,
                        ",".join(sorted(go_terms)),
                        ",".join(fragments),
                    )
                )

    return hits


def parse_jsonl(file: Path, extended: bool, repr_locs_only: bool) -> set[Hit]:
    hits = set()
    with file.open("rt") as fh:
        for line in map(str.rstrip, fh):
            data = json.loads(line)
            hits |= parse_external_json(
                data["results"], extended=extended, repr_locs_only=repr_locs_only
            )

    return hits


def format_xml(file: Path, extended: bool) -> set[Hit]:
    hits = set()

    for _, protein in ET.iterparse(file, events=("end",)):
        if protein.tag.rsplit("}", maxsplit=1)[-1] != "protein":
            continue

        seq_ids = []
        matches = None
        for child in protein:
            tag = child.tag.rsplit("}", maxsplit=1)[-1]
            if tag == "xref":
                seq_id = child.get("id")
                if seq_id is not None:
                    seq_ids.append(seq_id)
            elif tag == "matches":
                matches = child

        if not seq_ids or matches is None:
            protein.clear()
            continue

        for match in matches:
            if match.tag.rsplit("}", maxsplit=1)[-1] != "match":
                continue

            source = match.get("source")
            assert source is not None

            sig_acc = None
            locations = None
            for child in match:
                tag = child.tag.rsplit("}", maxsplit=1)[-1]
                if tag == "signature":
                    sig_acc = child.get("ac")
                elif tag == "locations":
                    locations = child

            if sig_acc is None or locations is None:
                continue

            for loc in locations:
                if loc.tag.rsplit("}", maxsplit=1)[-1] != "location":
                    continue

                fragments = []
                if extended:
                    for child in loc:
                        if (
                            child.tag.rsplit("}", maxsplit=1)[-1]
                            != "location-fragments"
                        ):
                            continue

                        for f in child:
                            if f.tag.rsplit("}", maxsplit=1)[-1] != "fragment":
                                continue

                            fragments.append(
                                (
                                    int(f.attrib["start"]),
                                    int(f.attrib["end"]),
                                    f.attrib["dc-status"],
                                )
                            )

                    fragments.sort(key=lambda x: (x[0], x[1]))
                    region = ",".join(
                        f"{start}-{end}-{dc}" for start, end, dc in fragments
                    )
                else:
                    start = int(loc.attrib["start"])
                    end = int(loc.attrib["end"])
                    region = f"{start}-{end}"

                for seq_id in seq_ids:
                    # TODO
                    hits.add((seq_id, source, sig_acc, region))

        protein.clear()

    return hits


if __name__ == "__main__":
    main()
