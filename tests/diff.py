#!/usr/bin/env python

import argparse
import gzip
import json
from pathlib import Path
import shutil
from tempfile import TemporaryDirectory
from typing import Any
import xml.etree.ElementTree as ET


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("expected", type=Path)
    parser.add_argument("observed", type=Path)

    args = parser.parse_args()

    with TemporaryDirectory() as tmpdir:
        exp_file = decompress_if_gzipped(args.expected, tmpdir)
        obs_file = decompress_if_gzipped(args.observed, tmpdir)

        extended = supports_extended(exp_file) and supports_extended(obs_file)

        exp_hits = parse(exp_file, extended=extended)
        obs_hits = parse(obs_file, extended=extended)

    assert exp_hits == obs_hits


def supports_extended(file: Path) -> bool:
    return file.name.endswith((".json", ".jsonl", ".xml"))


def decompress_if_gzipped(file: Path, tmpdir: TemporaryDirectory) -> Path:
    if file.suffix != ".gz":
        return file

    tmpfile = Path(tmpdir.name) / file.with_suffix("").name
    with gzip.open(file, "rb") as src, tmpfile.open("wb") as dst:
        shutil.copyfileobj(src, dst, length=1024 * 1024)

    return tmpfile


def parse(file: Path, extended: bool) -> set[tuple[str, str, str, str]]:
    fmt = file.suffix[1:]
    if fmt == "gff3":
        return parse_gff3(file)
    elif fmt == "json":
        return parse_json(file, extended=extended)
    elif fmt == "jsonl":
        return parse_jsonl(file, extended=extended)
    elif fmt == "tsv":
        return format_tsv(file)
    elif fmt == "xml":
        return format_xml(file, extended=extended)
    else:
        raise NotImplementedError(f"Format not supported: {fmt}")


def parse_gff3(file: Path) -> set[tuple[str, str, str, str]]:
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

            d_attrs = {}
            for attr in attrs.split(";"):
                key, value = attr.split("=")
                d_attrs[key] = value

            if source == "InterPro-N":
                _, sig_acc = d_attrs["Dbxref"].split(":", maxsplit=1)
            else:
                sig_acc = d_attrs.get("Alias", d_attrs["Name"])

            hits.add((seq_id, source, sig_acc, f"{start}-{end}"))

    return hits


def format_tsv(file: str) -> set[tuple[str, str, str, str]]:
    hits = set()
    with file.open("rt") as fh:
        for line in map(str.rstrip, fh):
            values = line.split("\t")
            assert len(values) == 15
            seq_id = values[0]
            sig_acc = values[4]
            start = int(values[6])
            end = values[7]
            source = values[9]
            hits.add((seq_id, source, sig_acc, f"{start}-{end}"))

    return hits


def parse_json(file: Path, extended: bool) -> set[tuple[str, str, str, str]]:
    with file.open("rt") as fh:
        data = json.load(fh)

    if "results" in data:
        return parse_external_json(data["results"], extended=extended)
    else:
        return parse_internal_json(data, extended=extended)


def parse_external_json(
    data: list[dict[str, Any]], extended: bool
) -> set[tuple[str, str, str, str]]:
    hits = set()
    for seq in data:
        for xref in seq["xref"]:
            seq_id = xref["id"]
            for match in seq["matches"]:
                source = match["source"]
                sig_acc = match["signature"]["accession"]
                for loc in match["locations"]:
                    fragments = []
                    if extended:
                        for f in sorted(
                            loc["location-fragments"], key=lambda x: (x["start"], x["end"])
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

                    hits.add((seq_id, source, sig_acc, ",".join(fragments)))

    return hits


def parse_internal_json(
    data: dict[str, dict[str, dict[str, Any]]], extended: bool
) -> set[tuple[str, str, str, str]]:
    hits = set()
    for seq_id, seq_matches in data.items():
        if not seq_matches:
            # No matches found in this sequence
            hits.add((seq_id, "N/A", "N/A", ""))
            continue

        for match in seq_matches.values():
            for loc in match["locations"]:
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
                        ",".join(fragments),
                    )
                )

    return hits


def parse_jsonl(file: Path, extended: bool) -> set[tuple[str, str, str, str]]:
    hits = set()
    with file.open("rt") as fh:
        for line in map(str.rstrip, fh):
            data = json.loads(line)
            hits |= parse_external_json(data["results"], extended=extended)

    return hits


def format_xml(file: Path, extended: bool) -> set[tuple[str, str, str, str]]:
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
                        if child.tag.rsplit("}", maxsplit=1)[-1] != "location-fragments":
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
                    region = ",".join(f"{start}-{end}-{dc}" for start, end, dc in fragments)
                else:
                    start = int(loc.attrib["start"])
                    end = int(loc.attrib["end"])
                    region = f"{start}-{end}"

                for seq_id in seq_ids:
                    hits.add((seq_id, source, sig_acc, region))

        protein.clear()

    return hits


if __name__ == "__main__":
    main()
