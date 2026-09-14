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

    spec = importlib.util.spec_from_file_location(
        "editdistance_cigar_filter",
        SCRIPT_PATH,
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


HEADER = [
    "read_id", "RNAseq_SE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags",
]


ROWS = [
    # Baseline: valid RNAseq_SE read with no splice gap.
    [
        "good_read", "U",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # NM failure.
    [
        "high_nm", "U",
        "chr1", 100, 119, "+", "20M", 99, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # MAPQ failure.
    [
        "low_mapq", "U",
        "chr1", 100, 119, "+", "20M", 0, 0,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # One N is allowed. Left aligned block is longer, so final coordinates
    # should be trimmed to the left block.
    [
        "one_N_left", "U",
        "chr1", 100, 227, "+", "20M100N8M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # One N is allowed. Right aligned block is longer, so final coordinates
    # should be trimmed to the right block.
    [
        "one_N_right", "U",
        "chr1", 100, 227, "+", "8M100N20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # More than one N must be rejected.
    [
        "multi_N", "U",
        "chr1", 100, 134, "+", "5M10N5M10N5M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # In NM + clipping mode, 10 clipped bases must exceed threshold=2.
    [
        "softclip_penalty", "U",
        "chr1", 100, 109, "+", "10S10M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Small clipping still fits threshold=2 and should pass.
    [
        "small_softclip", "U",
        "chr1", 100, 117, "+", "2S18M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],
]


EXPECTED_FILTERED = {
    "good_read",
    "one_N_left",
    "one_N_right",
    "small_softclip",
}

EXPECTED_OUT = {
    "high_nm",
    "low_mapq",
    "multi_N",
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

    with tempfile.TemporaryDirectory(prefix="rnaseq_se_filtering_") as tmp:
        tmp_path = Path(tmp)
        input_name = "rnaseq_se_test.tab"
        input_path = tmp_path / input_name
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        write_tsv(input_path, HEADER, ROWS)

        module.editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            2,                      # r1 edit-distance threshold
            2,                      # unused for RNAseq_SE
            10,                     # r1 MAPQ threshold
            10,                     # unused for RNAseq_SE
            200,                    # unused for RNAseq_SE
            "no",
            "not explorer",
            "RNAseq_SE",
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

        # Unified output format: RNA half populated, DNA half empty.
        baseline = read_row_by_id(filtered_path, "good_read")
        if baseline is None:
            failures.append("Cannot find good_read in filtered output.")
        else:
            if baseline.get("rna_chr") != "chr1":
                failures.append(
                    f"Unexpected RNA chromosome: {baseline.get('rna_chr')}"
                )
            if baseline.get("rna_start") != "100":
                failures.append(
                    f"Unexpected RNA start: {baseline.get('rna_start')}"
                )
            if baseline.get("rna_end") != "119":
                failures.append(
                    f"Unexpected RNA end: {baseline.get('rna_end')}"
                )

            dna_columns = [
                "dna_chr", "dna_start", "dna_end", "dna_strand",
                "dna_cigar", "dna_NM", "dna_mapq",
                "dna_secondary_alignments", "dna_other_tags",
            ]
            non_star_dna = {
                column: baseline.get(column)
                for column in dna_columns
                if baseline.get(column) != "*"
            }
            if non_star_dna:
                failures.append(
                    "RNAseq_SE filtered output unexpectedly contains DNA fields: "
                    f"{non_star_dna}"
                )

        # One-N trimming: longer left block should be selected.
        left = read_row_by_id(filtered_path, "one_N_left")
        if left is None:
            failures.append("Cannot find one_N_left in filtered output.")
        else:
            if left.get("rna_start") != "100" or left.get("rna_end") != "119":
                failures.append(
                    "one_N_left coordinates mismatch: "
                    f"{left.get('rna_start')}-{left.get('rna_end')} "
                    "(expected 100-119)"
                )

        # One-N trimming: longer right block should be selected.
        right = read_row_by_id(filtered_path, "one_N_right")
        if right is None:
            failures.append("Cannot find one_N_right in filtered output.")
        else:
            if right.get("rna_start") != "208" or right.get("rna_end") != "227":
                failures.append(
                    "one_N_right coordinates mismatch: "
                    f"{right.get('rna_start')}-{right.get('rna_end')} "
                    "(expected 208-227)"
                )

        if failures:
            print("FAIL: RNAseq_SE integration test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: RNAseq_SE integration test")
        print(f"  filtered: {len(filtered_ids)}/{len(ROWS)}")
        print(f"  out:      {len(out_ids)}/{len(ROWS)}")
        print("  high NM:             rejected as expected")
        print("  low MAPQ:            rejected as expected")
        print("  one N:               accepted as expected")
        print("  one-N coordinates:   trimmed to the longer aligned block")
        print("  multiple N:          rejected as expected")
        print("  clipping penalty:    applied in NM + N_softClipp_bp mode")
        print("  output format:       RNA half populated, DNA half empty")


if __name__ == "__main__":
    main()
