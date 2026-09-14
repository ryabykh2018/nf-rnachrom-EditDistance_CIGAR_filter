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
    "read_id", "OTA_PE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags",
]

ROWS = [
    ["good_pair", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 150, 169, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_high_nm", "UU",
     "chr1", 100, 119, "+", "20M", 99, 60,
     "chr1", 150, 169, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r2_high_nm", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 150, 169, "-", "20M", 99, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_low_mapq", "UU",
     "chr1", 100, 119, "+", "20M", 0, 0,
     "chr1", 150, 169, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r2_low_mapq", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 150, 169, "-", "20M", 0, 0,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r1_has_N", "UU",
     "chr1", 100, 219, "+", "10M100N10M", 0, 60,
     "chr1", 250, 269, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["r2_has_N", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 150, 269, "-", "10M100N10M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["different_chr", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr2", 150, 169, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    ["too_far", "UU",
     "chr1", 100, 119, "+", "20M", 0, 60,
     "chr1", 1000, 1019, "-", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],

    # Same strand is intentional: strand orientation is an upstream
    # proper-pair requirement, not a downstream OTA_PE filter here.
    ["nearby_pair", "UU",
     "chr1", 500, 519, "+", "20M", 0, 60,
     "chr1", 650, 669, "+", "20M", 0, 60,
     "*", "*", "NH:i:1", "NH:i:1"],
]

EXPECTED_FILTERED = {"good_pair", "nearby_pair"}
EXPECTED_OUT = {
    "r1_high_nm", "r2_high_nm",
    "r1_low_mapq", "r2_low_mapq",
    "r1_has_N", "r2_has_N",
    "different_chr", "too_far",
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

    with tempfile.TemporaryDirectory(prefix="ota_pe_filtering_") as tmp:
        tmp_path = Path(tmp)
        input_name = "ota_pe_test.tab"
        input_path = tmp_path / input_name
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        write_tsv(input_path, HEADER, ROWS)

        module.editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            2, 2,
            10, 10,
            200,
            "no",
            "not explorer",
            "OTA_PE",
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

        baseline = read_row_by_id(filtered_path, "good_pair")
        if baseline is None:
            failures.append("Cannot find good_pair in filtered output.")
        else:
            if baseline.get("rna_chr") != "*":
                failures.append("OTA_PE output should leave the RNA half empty.")
            if baseline.get("dna_chr") != "chr1":
                failures.append(
                    f"Unexpected merged DNA chromosome: {baseline.get('dna_chr')}"
                )
            if baseline.get("dna_start") != "100":
                failures.append(
                    f"Unexpected merged DNA start: {baseline.get('dna_start')}"
                )
            if baseline.get("dna_end") != "169":
                failures.append(
                    f"Unexpected merged DNA end: {baseline.get('dna_end')}"
                )

        if failures:
            print("FAIL: OTA_PE integration test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: OTA_PE integration test")
        print(f"  filtered: {len(filtered_ids)}/{len(ROWS)}")
        print(f"  out:      {len(out_ids)}/{len(ROWS)}")
        print("  r1/r2 NM failures:        rejected as expected")
        print("  r1/r2 MAPQ failures:      rejected as expected")
        print("  r1/r2 N-containing CIGAR: rejected as expected")
        print("  different chromosome:     rejected as expected")
        print("  excessive distance:       rejected as expected")
        print("  accepted mates:           merged into one DNA interval")
        print("  strand orientation:       intentionally not tested downstream")


if __name__ == "__main__":
    main()
