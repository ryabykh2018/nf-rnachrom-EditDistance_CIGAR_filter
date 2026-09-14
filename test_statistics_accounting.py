#!/usr/bin/env python3

import os
import sys
import tempfile
import importlib.util
from pathlib import Path

import pandas as pd

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
    # Passes all filters.
    [
        "good_read", "U",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Normal filtering reject: edit distance exceeds threshold.
    [
        "high_nm", "U",
        "chr1", 100, 119, "+", "20M", 99, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Normal filtering reject: valid MAPQ, but below threshold.
    [
        "low_mapq", "U",
        "chr1", 100, 119, "+", "20M", 0, 5,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Normal filtering reject: RNAseq_SE allows at most one N.
    [
        "multi_N", "U",
        "chr1", 100, 134, "+", "5M10N5M10N5M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Validation reject: unsupported CIGAR operation.
    [
        "invalid_cigar", "U",
        "chr1", 100, 119, "+", "20Q", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Validation reject: NM is not numeric.
    [
        "invalid_nm", "U",
        "chr1", 100, 119, "+", "20M", "abc", 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Validation reject: MAPQ parser intentionally rejects "60.0".
    [
        "invalid_mapq", "U",
        "chr1", 100, 119, "+", "20M", 0, "60.0",
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],

    # Validation reject: genomic coordinates must be positive.
    [
        "invalid_coordinates", "U",
        "chr1", 0, 19, "+", "20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "NH:i:1", "*",
    ],
]


EXPECTED_FILTERED_IDS = {"good_read"}

EXPECTED_NORMAL_REJECT_IDS = {
    "high_nm",
    "low_mapq",
    "multi_N",
}

EXPECTED_VALIDATION_REJECT_IDS = {
    "invalid_cigar",
    "invalid_nm",
    "invalid_mapq",
    "invalid_coordinates",
}

EXPECTED_VALIDATION_COUNTS = {
    "invalid_CIGAR_r1": 1,
    "invalid_NM": 1,
    "invalid_MAPQ": 1,
    "invalid_coordinates": 1,
}


def write_tsv(path, header, rows):
    """Write a small tab-separated test input."""
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def read_ids(path):
    """Return data-row IDs from a tab-separated output file."""
    ids = set()
    with open(path) as f:
        next(f, None)
        for line in f:
            if line.strip():
                ids.add(line.split("\t", 1)[0])
    return ids


def count_data_rows(path):
    """Count non-header data rows."""
    with open(path) as f:
        next(f, None)
        return sum(1 for line in f if line.strip())


def read_stat_table(path):
    """Read a statistics table, preserving an empty table with headers."""
    return pd.read_csv(path, sep="\t")


def stat_sum(df):
    """Return the sum of the N column."""
    if df.empty:
        return 0
    return int(df["N"].sum())


def main():
    module = load_filter_module()
    failures = []

    with tempfile.TemporaryDirectory(prefix="statistics_accounting_") as tmp:
        tmp_path = Path(tmp)
        input_name = "statistics_test.tab"
        input_path = tmp_path / input_name
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        write_tsv(input_path, HEADER, ROWS)

        module.editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            2,      # r1 edit-distance threshold
            2,      # unused for RNAseq_SE
            10,     # r1 MAPQ threshold
            10,     # unused for RNAseq_SE
            200,    # unused for RNAseq_SE
            "no",
            "not explorer",
            "RNAseq_SE",
            input_name,
            str(tmp_path) + os.sep,
            str(output_dir) + os.sep,
        )

        filtered_path = output_dir / f"filtered_{input_name}"
        out_path = output_dir / f"out_{input_name}"
        cigar_filtered_path = output_dir / f"cigar_stat_filtered_{input_name}"
        cigar_out_path = output_dir / f"cigar_stat_out_{input_name}"
        validation_path = output_dir / f"validation_reject_stat_out_{input_name}"

        required_files = [
            filtered_path,
            out_path,
            cigar_filtered_path,
            cigar_out_path,
            validation_path,
        ]

        for path in required_files:
            if not path.exists():
                failures.append(f"Missing expected output file: {path.name}")

        if failures:
            print("FAIL: statistics accounting test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        filtered_ids = read_ids(filtered_path)
        out_ids = read_ids(out_path)

        expected_out_ids = (
            EXPECTED_NORMAL_REJECT_IDS
            | EXPECTED_VALIDATION_REJECT_IDS
        )

        if filtered_ids != EXPECTED_FILTERED_IDS:
            failures.append(
                "Filtered IDs mismatch\n"
                f"  expected: {sorted(EXPECTED_FILTERED_IDS)}\n"
                f"  actual:   {sorted(filtered_ids)}"
            )

        if out_ids != expected_out_ids:
            failures.append(
                "Out IDs mismatch\n"
                f"  expected: {sorted(expected_out_ids)}\n"
                f"  actual:   {sorted(out_ids)}"
            )

        cigar_filtered = read_stat_table(cigar_filtered_path)
        cigar_out = read_stat_table(cigar_out_path)
        validation = read_stat_table(validation_path)

        filtered_rows = count_data_rows(filtered_path)
        out_rows = count_data_rows(out_path)

        filtered_cigar_n = stat_sum(cigar_filtered)
        out_cigar_n = stat_sum(cigar_out)
        validation_n = stat_sum(validation)

        # Every filtered contact must be represented exactly once
        # in filtered CIGAR statistics.
        if filtered_rows != filtered_cigar_n:
            failures.append(
                "Filtered accounting mismatch\n"
                f"  filtered rows:        {filtered_rows}\n"
                f"  filtered CIGAR stats: {filtered_cigar_n}"
            )

        # Every rejected contact must be represented exactly once either
        # as a normal CIGAR-based reject or as a validation reject.
        if out_rows != out_cigar_n + validation_n:
            failures.append(
                "Out accounting mismatch\n"
                f"  out rows:             {out_rows}\n"
                f"  out CIGAR stats:      {out_cigar_n}\n"
                f"  validation rejects:   {validation_n}\n"
                f"  combined stats total: {out_cigar_n + validation_n}"
            )

        # This synthetic dataset has exactly three normal filtering rejects.
        if out_cigar_n != len(EXPECTED_NORMAL_REJECT_IDS):
            failures.append(
                "Unexpected number of normal filtering rejects in "
                "cigar_stat_out\n"
                f"  expected: {len(EXPECTED_NORMAL_REJECT_IDS)}\n"
                f"  actual:   {out_cigar_n}"
            )

        # Validation categories must not leak into cigar_stat_out.
        cigar_type_column = next(
            (col for col in cigar_out.columns if col.startswith("CIGAR type")),
            None,
        )

        if cigar_type_column is None:
            failures.append(
                "Cannot find the CIGAR type column in cigar_stat_out."
            )
        else:
            leaked = [
                str(value)
                for value in cigar_out[cigar_type_column].tolist()
                if str(value).startswith("invalid")
            ]
            if leaked:
                failures.append(
                    "Validation categories leaked into cigar_stat_out: "
                    + ", ".join(leaked)
                )

        # Validation file must contain only the expected technical rejects.
        actual_validation_counts = {}
        if not validation.empty:
            actual_validation_counts = dict(
                zip(
                    validation["validation_reject_reason"],
                    validation["N"].astype(int),
                )
            )

        if actual_validation_counts != EXPECTED_VALIDATION_COUNTS:
            failures.append(
                "Validation rejection counts mismatch\n"
                f"  expected: {EXPECTED_VALIDATION_COUNTS}\n"
                f"  actual:   {actual_validation_counts}"
            )

        # The validation and normal-reject counts are deliberately disjoint.
        if validation_n != len(EXPECTED_VALIDATION_REJECT_IDS):
            failures.append(
                "Unexpected number of validation rejects\n"
                f"  expected: {len(EXPECTED_VALIDATION_REJECT_IDS)}\n"
                f"  actual:   {validation_n}"
            )

        if failures:
            print("FAIL: statistics accounting test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: statistics accounting test")
        print(f"  filtered rows:              {filtered_rows}")
        print(f"  filtered CIGAR stats:       {filtered_cigar_n}")
        print(f"  out rows:                   {out_rows}")
        print(f"  normal CIGAR rejects:       {out_cigar_n}")
        print(f"  validation rejects:         {validation_n}")
        print(
            "  out accounting:             "
            f"{out_cigar_n} + {validation_n} = {out_rows}"
        )
        print("  invalid_* in cigar_stat_out: absent")
        print("  validation reasons:         separated correctly")


if __name__ == "__main__":
    main()
