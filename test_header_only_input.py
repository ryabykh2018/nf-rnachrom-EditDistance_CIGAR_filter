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


def count_data_rows(path):
    """Count non-header, non-empty rows."""
    with open(path) as f:
        next(f, None)
        return sum(1 for line in f if line.strip())


def main():
    module = load_filter_module()
    failures = []

    with tempfile.TemporaryDirectory(prefix="header_only_") as tmp:
        tmp_path = Path(tmp)
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        input_name = "header_only.tab"
        input_path = tmp_path / input_name

        # Valid input table with a header and zero contacts.
        input_path.write_text("\t".join(HEADER) + "\n")

        try:
            module.editDistance_and_CIGAR_filter(
                "NM + N_softClipp_bp",
                2,
                2,
                10,
                10,
                200,
                "no",
                "not explorer",
                "RNAseq_SE",
                input_name,
                str(tmp_path) + os.sep,
                str(output_dir) + os.sep,
            )
        except Exception as exc:
            print("FAIL: header-only input test")
            print(
                "  editDistance_and_CIGAR_filter raised an exception for "
                f"a valid header-only input: {type(exc).__name__}: {exc}"
            )
            sys.exit(1)

        filtered_path = output_dir / f"filtered_{input_name}"
        out_path = output_dir / f"out_{input_name}"
        cigar_filtered_path = output_dir / f"cigar_stat_filtered_{input_name}"
        cigar_out_path = output_dir / f"cigar_stat_out_{input_name}"
        validation_path = (
            output_dir / f"validation_reject_stat_out_{input_name}"
        )

        required_files = [
            filtered_path,
            out_path,
            cigar_filtered_path,
            cigar_out_path,
            validation_path,
        ]

        for required_path in required_files:
            if not required_path.exists():
                failures.append(
                    f"Missing expected output file: {required_path.name}"
                )

        if filtered_path.exists():
            if count_data_rows(filtered_path) != 0:
                failures.append(
                    "filtered_* contains data rows for a header-only input."
                )

            lines = filtered_path.read_text().splitlines()
            if len(lines) != 1:
                failures.append(
                    "filtered_* should contain exactly one header line."
                )

        if out_path.exists():
            if count_data_rows(out_path) != 0:
                failures.append(
                    "out_* contains data rows for a header-only input."
                )

            lines = out_path.read_text().splitlines()
            if len(lines) != 1:
                failures.append(
                    "out_* should contain exactly one header line."
                )

        for stat_path in [
            cigar_filtered_path,
            cigar_out_path,
            validation_path,
        ]:
            if not stat_path.exists():
                continue

            try:
                df = pd.read_csv(stat_path, sep="\t")
            except Exception as exc:
                failures.append(
                    f"{stat_path.name} is not a readable TSV statistics file: "
                    f"{type(exc).__name__}: {exc}"
                )
                continue

            if "N" not in df.columns:
                failures.append(
                    f"{stat_path.name} does not contain the expected N column."
                )
            elif not df.empty and int(df["N"].sum()) != 0:
                failures.append(
                    f"{stat_path.name} contains non-zero statistics for "
                    "a header-only input."
                )

        uca_path = output_dir / f"id_reads_for_ucaRNAs_{input_name}"
        if uca_path.exists():
            failures.append(
                "RNAseq_SE unexpectedly created id_reads_for_ucaRNAs_*."
            )

        if failures:
            print("FAIL: header-only input test")
            for failure in failures:
                print("\n" + failure)
            sys.exit(1)

        print("PASS: header-only input test")
        print("  input:                    valid header, zero contacts")
        print("  filtered data rows:       0")
        print("  out data rows:            0")
        print("  CIGAR statistics counts:  0")
        print("  validation reject counts: 0")
        print("  required output files:    created")
        print("  ucaRNA output:            absent for RNAseq_SE")


if __name__ == "__main__":
    main()
