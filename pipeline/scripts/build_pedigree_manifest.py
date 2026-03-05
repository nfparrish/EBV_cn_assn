#!/usr/bin/env python3
"""Build a pedigree-annotated manifest from a 1000G .ped file and CRAM index.

Usage:
    python build_pedigree_manifest.py <ped_file> <output_csv> [options]

Options:
    --cram-base-url <url>   Base URL for 30x CRAMs (default: FTP URL for 1000G high-coverage)
    --cram-local-dir <dir>  Local directory with CRAMs (overrides --cram-base-url)
    --max-trios <n>         Maximum number of complete trios to include (default: all)
    --populations <pop,...>  Comma-separated population filter (default: all)
    --require-complete-trio  Only include samples that are part of a complete trio
"""
import argparse
import csv
import pathlib
import sys

DEFAULT_CRAM_BASE = (
    "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/data"
)
DEFAULT_CRAM_SUFFIX = ".final.cram"


def parse_ped(ped_path):
    """Parse a 1000 Genomes .ped file.
    Returns list of dicts with family_id, sample_id, father_id, mother_id, sex, population, relationship.
    """
    records = []
    with open(ped_path, "r") as fh:
        header_line = fh.readline()
        for line in fh:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            family_id = parts[0]
            sample_id = parts[1]
            father_id = parts[2]
            mother_id = parts[3]
            sex_code = parts[4]
            population = parts[6]
            relationship = parts[7] if len(parts) > 7 else ""

            sex = "M" if sex_code == "1" else ("F" if sex_code == "2" else "U")
            records.append({
                "family_id": family_id,
                "sample_id": sample_id,
                "father_id": father_id if father_id != "0" else "",
                "mother_id": mother_id if mother_id != "0" else "",
                "sex": sex,
                "population": population,
                "relationship": relationship.strip().lower() if relationship else "",
            })
    return records


def build_cram_path(sample_id, population, base_url, local_dir=None):
    """Construct the CRAM path for a given sample."""
    if local_dir:
        local = pathlib.Path(local_dir) / f"{sample_id}{DEFAULT_CRAM_SUFFIX}"
        return str(local)
    return f"{base_url}/{population}/{sample_id}/{sample_id}{DEFAULT_CRAM_SUFFIX}"


def identify_complete_trios(records):
    """Return list of (child, father, mother) tuples for complete trios."""
    by_id = {r["sample_id"]: r for r in records}
    trios = []
    for rec in records:
        if rec["father_id"] and rec["mother_id"]:
            if rec["father_id"] in by_id and rec["mother_id"] in by_id:
                trios.append((rec, by_id[rec["father_id"]], by_id[rec["mother_id"]]))
    return trios


def main():
    parser = argparse.ArgumentParser(description="Build pedigree-annotated manifest")
    parser.add_argument("ped_file", help="Path to 1000G .ped file")
    parser.add_argument("output_csv", help="Output manifest CSV path")
    parser.add_argument("--cram-base-url", default=DEFAULT_CRAM_BASE)
    parser.add_argument("--cram-local-dir", default=None)
    parser.add_argument("--max-trios", type=int, default=0, help="Max trios (0=all)")
    parser.add_argument("--populations", default="", help="Comma-separated population filter")
    parser.add_argument("--require-complete-trio", action="store_true")
    args = parser.parse_args()

    records = parse_ped(args.ped_file)
    pop_filter = set(args.populations.split(",")) if args.populations else set()

    if pop_filter:
        records = [r for r in records if r["population"] in pop_filter]

    trios = identify_complete_trios(records)

    if args.max_trios and args.max_trios > 0:
        trios = trios[: args.max_trios]

    if args.require_complete_trio:
        trio_sample_ids = set()
        for child, father, mother in trios:
            trio_sample_ids.update([child["sample_id"], father["sample_id"], mother["sample_id"]])
        records = [r for r in records if r["sample_id"] in trio_sample_ids]
    elif args.max_trios and args.max_trios > 0:
        trio_sample_ids = set()
        for child, father, mother in trios:
            trio_sample_ids.update([child["sample_id"], father["sample_id"], mother["sample_id"]])
        records = [r for r in records if r["sample_id"] in trio_sample_ids]

    # Deduplicate (parents can appear in multiple trios)
    seen = set()
    unique_records = []
    for rec in records:
        if rec["sample_id"] not in seen:
            seen.add(rec["sample_id"])
            unique_records.append(rec)
    records = unique_records

    header = [
        "sample_id", "dataset", "population", "sex", "sequence_lab",
        "cram_path", "family_id", "father_id", "mother_id", "relationship",
    ]

    with open(args.output_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=header)
        writer.writeheader()
        for rec in sorted(records, key=lambda r: (r["family_id"], r["sample_id"])):
            cram_path = build_cram_path(
                rec["sample_id"], rec["population"],
                args.cram_base_url, args.cram_local_dir,
            )
            writer.writerow({
                "sample_id": rec["sample_id"],
                "dataset": "1000G_PHASE3",
                "population": rec["population"],
                "sex": rec["sex"],
                "sequence_lab": "NYGC",
                "cram_path": cram_path,
                "family_id": rec["family_id"],
                "father_id": rec["father_id"] or "0",
                "mother_id": rec["mother_id"] or "0",
                "relationship": rec["relationship"] or "unrelated",
            })

    n_trios = len(trios)
    n_samples = len(records)
    print(f"Wrote {n_samples} samples ({n_trios} complete trios) to {args.output_csv}")


if __name__ == "__main__":
    main()
