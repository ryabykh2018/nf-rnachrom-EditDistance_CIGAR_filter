#!/usr/bin/env python3

import os
import sys
import tempfile
import importlib.util
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")
SCRIPT_PATH = Path(__file__).with_name("EditDistance_CIGAR_filter.py")


def load_filter_module():
    """Load EditDistance_CIGAR_filter.py from the same directory as this test."""
    if not SCRIPT_PATH.exists():
        raise FileNotFoundError(
            f"Cannot find {SCRIPT_PATH.name} next to {Path(__file__).name}"
        )
    spec = importlib.util.spec_from_file_location("editdistance_cigar_filter", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HEADER = [
    "read_id", "OTA_SE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags",
]

ROWS = [
    # Baseline: valid OTA_SE read.
    ["good_read", "U",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],

    # NM failure.
    ["high_nm", "U",
     "chr1", 100, 119, "+", "20M", 99, 60,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],

    # MAPQ failure.
    ["low_mapq", "U",
     "chr1", 100, 119, "+", "20M", 0, 0,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],

    # OTA_SE forbids any splice N.
    ["has_N", "U",
     "chr1", 100, 219, "+", "10M100N10M", 0, 60,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],

    # With NM + clipping mode, 10 clipped bases should exceed threshold=2.
    ["softclip_penalty", "U",
     "chr1", 100, 109, "+", "10S10M", 0, 60,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],

    # Small clipping that still fits threshold=2 should pass.
    ["small_softclip", "U",
     "chr1", 100, 117, "+", "2S18M", 0, 60,
     "*", "*", "*", "*", "*", "*", "*",
     "*", "*", "NH:i:1", "*"],
]

EXPECTED_FILTERED = {
    "good_read",
    "small_softclip",
}

EXPECTED_OUT = {
    "high_nm",
    "low_mapq",
    "has_N",
    "softclip_penalty",
}


def write_tsv(path, header, rows):
    """Write a small tab-separated test input."""
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def read_ids(path):
    """Return all read IDs from an output table, excluding the header."""
    ids = set()
    with open(path) as f:
        next(f, None)
        for line in f:
            if line.strip():
                ids.add(line.split("\t", 1)[0])
    return ids


def read_row_by_id(path, read_id):
    """Read one output row and return it as a dictionary."""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line in f:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == read_id:
                return dict(zip(header, fields))
    return None


def main():
    module = load_filter_module()

    with tempfile.TemporaryDirectory(prefix="ota_se_filtering_") as tmp:
        tmp_path = Path(tmp)
        input_name = "ota_se_test.tab"
        input_path = tmp_path / input_name
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        write_tsv(input_path, HEADER, ROWS)

        module.editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            2,                      # r1 edit-distance threshold
            2,                      # unused for OTA_SE
            10,                     # r1 MAPQ threshold
            10,                     # unused for OTA_SE
            200,                    # unused for OTA_SE
            "no",
            "not explorer",
            "OTA_SE",
            input_name,
            str(tmp_path) + os.sep,
            str(output_dir) + os.sep,
        )

        filtered_path = output_dir / f"filtered_{input_name}"
        out_path = output_dir / f"out_{input_name}"

        filtered_ids = read_ids(filtered_path)
        out_ids = read_ids(out_path)

        failures = []

        if filtered_ids != EXPECTED_FILTERED:
            failures.append(
                "filtered IDs mismatch\n"
                f"  expected: {sorted(EXPECTED_FILTERED)}\n"
                f"  actual:   {sorted(filtered_ids)}"
            )

        if out_ids != EXPECTED_OUT:
            failures.append(
                "out IDs mismatch\n"
                f"  expected: {sorted(EXPECTED_OUT)}\n"
                f"  actual:   {sorted(out_ids)}"
            )

        baseline = read_row_by_id(filtered_path, "good_read")
        if baseline is None:
            failures.append("Cannot find good_read in filtered output.")
        else:
            # Unified output format: RNA half must be empty for OTA_SE.
            rna_columns = [
                "rna_chr", "rna_start", "rna_end", "rna_strand",
                "rna_cigar", "rna_NM", "rna_mapq",
                "rna_secondary_alignments", "rna_other_tags",
            ]
            non_star_rna = {
                column: baseline.get(column)
                for column in rna_columns
                if baseline.get(column) != "*"
            }
            if non_star_rna:
                failures.append(
                    "OTA_SE filtered output unexpectedly contains RNA fields: "
                    f"{non_star_rna}"
                )

            if baseline.get("dna_chr") != "chr1":
                failures.append(
                    f"Unexpected DNA chromosome: {baseline.get('dna_chr')}"
                )
            if baseline.get("dna_start") != "100":
                failures.append(
                    f"Unexpected DNA start: {baseline.get('dna_start')}"
                )
            if baseline.get("dna_end") != "119":
                failures.append(
                    f"Unexpected DNA end: {baseline.get('dna_end')}"
                )

        if failures:
            print("FAIL: OTA_SE integration test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: OTA_SE integration test")
        print(f"  filtered: {len(filtered_ids)}/{len(ROWS)}")
        print(f"  out:      {len(out_ids)}/{len(ROWS)}")
        print("  high NM:            rejected as expected")
        print("  low MAPQ:           rejected as expected")
        print("  N-containing CIGAR: rejected as expected")
        print("  clipping penalty:   applied in NM + N_softClipp_bp mode")
        print("  output format:      RNA half empty, DNA half populated")


if __name__ == "__main__":
    main()
