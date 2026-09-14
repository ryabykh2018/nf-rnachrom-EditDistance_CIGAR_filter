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
    "read_id", "RNAseq_PE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags",
]

ROWS = [
    ["r1_good_r2_good", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 200, 219, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_good_r2_high_nm", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 200, 219, "-", "20M", 99, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_good_r2_low_mapq", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 200, 219, "-", "20M", 0, 0,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_good_r2_multi_N", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 200, 234, "-", "5M10N5M10N5M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_good_r2_heavy_clip", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 200, 209, "-", "10S10M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_multi_N_r2_good", "UU",
     "chr1", 100, 134, "+", "5M10N5M10N5M", 0, 60,
     "chr1", 200, 219, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_high_nm_r2_good", "UU",
     "chr1", 100, 119, "+", "20M", 99, 60,
     "chr1", 200, 219, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_low_mapq_r2_good", "UU",
     "chr1", 100, 119, "+", "20M", 0, 0,
     "chr1", 200, 219, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],
]

EXPECTED_FILTERED = {
    "r1_good_r2_good",
    "r1_good_r2_high_nm",
    "r1_good_r2_low_mapq",
    "r1_good_r2_multi_N",
    "r1_good_r2_heavy_clip",
}

EXPECTED_OUT = {
    "r1_multi_N_r2_good",
    "r1_high_nm_r2_good",
    "r1_low_mapq_r2_good",
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

    with tempfile.TemporaryDirectory(prefix="rnaseq_pe_r1_only_") as tmp:
        tmp_path = Path(tmp)
        input_name = "rnaseq_pe_test.tab"
        input_path = tmp_path / input_name
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        write_tsv(input_path, HEADER, ROWS)

        module.editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            2,
            2,
            10,
            10,
            200,
            "no",
            "not explorer",
            "RNAseq_PE",
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

        row = read_row_by_id(filtered_path, "r1_good_r2_good")
        if row is None:
            failures.append("Cannot find baseline row in filtered output.")
        else:
            expected_r2_columns = [
                "dna_chr", "dna_start", "dna_end", "dna_strand",
                "dna_cigar", "dna_NM", "dna_mapq",
                "dna_secondary_alignments", "dna_other_tags",
            ]
            non_star = {
                column: row.get(column)
                for column in expected_r2_columns
                if row.get(column) != "*"
            }
            if non_star:
                failures.append(
                    "RNAseq_PE filtered output unexpectedly retains r2 fields: "
                    f"{non_star}"
                )

        if failures:
            print("FAIL: RNAseq_PE r1-only integration test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: RNAseq_PE r1-only integration test")
        print(f"  filtered: {len(filtered_ids)}/{len(ROWS)}")
        print(f"  out:      {len(out_ids)}/{len(ROWS)}")
        print("  r2 high NM:        ignored as expected")
        print("  r2 low MAPQ:       ignored as expected")
        print("  r2 multiple N:     ignored as expected")
        print("  r2 heavy clipping: ignored as expected")
        print("  r1 N/NM/MAPQ failures: rejected as expected")


if __name__ == "__main__":
    main()
