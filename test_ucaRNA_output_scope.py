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


ATA_HEADER = [
    "read_id", "ATA_pairtype",
    "rna_chr", "rna_start", "rna_end", "rna_strand",
    "rna_cigar", "rna_NM", "rna_mapq",
    "dna_chr", "dna_start", "dna_end", "dna_strand",
    "dna_cigar", "dna_NM", "dna_mapq",
    "rna_secondary_alignments", "dna_secondary_alignments",
    "rna_other_tags", "dna_other_tags",
]

ATA_ROW = [
    "ata_good", "UU",
    "chr1", 100, 119, "+", "20M", 0, 60,
    "chr1", 200, 219, "-", "20M", 0, 60,
    "*", "*", "NH:i:1", "NH:i:1",
]


RNASEQ_SE_HEADER = [
    "read_id", "RNAseq_SE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags",
]

RNASEQ_SE_ROW = [
    "rnaseq_se_good", "U",
    "chr1", 100, 119, "+", "20M", 0, 60,
    "*", "*", "*", "*", "*", "*", "*",
    "*", "*", "NH:i:1", "*",
]


RNASEQ_PE_HEADER = [
    "read_id", "RNAseq_PE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags",
]

RNASEQ_PE_ROW = [
    "rnaseq_pe_good", "UU",
    "chr1", 100, 119, "+", "20M", 0, 60,
    "chr1", 200, 219, "-", "20M", 0, 60,
    "*", "*", "NH:i:1", "NH:i:1",
]


OTA_SE_HEADER = [
    "read_id", "OTA_SE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags",
]

OTA_SE_ROW = [
    "ota_se_good", "U",
    "chr1", 100, 119, "+", "20M", 0, 60,
    "*", "*", "*", "*", "*", "*", "*",
    "*", "*", "NH:i:1", "*",
]


OTA_PE_HEADER = [
    "read_id", "OTA_PE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags",
]

OTA_PE_ROW = [
    "ota_pe_good", "UU",
    "chr1", 100, 119, "+", "20M", 0, 60,
    "chr1", 150, 169, "-", "20M", 0, 60,
    "*", "*", "NH:i:1", "NH:i:1",
]


def write_tsv(path, header, row):
    """Write a one-row tab-separated input file."""
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        f.write("\t".join(map(str, row)) + "\n")


def run_filter(module, tmp_path, experiment_type, header, row, assembly, suffix):
    """Run the filter once and return the expected ucaRNA output path."""
    input_name = f"{suffix}.tab"
    input_path = tmp_path / input_name
    output_dir = tmp_path / "output"
    output_dir.mkdir(exist_ok=True)

    write_tsv(input_path, header, row)

    module.editDistance_and_CIGAR_filter(
        "NM + N_softClipp_bp",
        2,
        2,
        10,
        10,
        200,
        assembly,
        "not explorer",
        experiment_type,
        input_name,
        str(tmp_path) + os.sep,
        str(output_dir) + os.sep,
    )

    return output_dir / f"id_reads_for_ucaRNAs_{input_name}"


def assert_file_exists(path, label, failures):
    if not path.exists():
        failures.append(
            f"{label}: expected id_reads_for_ucaRNAs_* to be created."
        )


def assert_file_absent(path, label, failures):
    if path.exists():
        failures.append(
            f"{label}: id_reads_for_ucaRNAs_* was created unexpectedly."
        )


def main():
    module = load_filter_module()
    failures = []

    with tempfile.TemporaryDirectory(prefix="ucarna_scope_") as tmp:
        tmp_path = Path(tmp)

        ata_yes = run_filter(
            module, tmp_path,
            "ATA, not iMARGI",
            ATA_HEADER, ATA_ROW,
            "yes",
            "ata_yes",
        )
        assert_file_exists(
            ata_yes,
            "ATA, not iMARGI + Assembly_of_ucaRNAs=yes",
            failures,
        )

        rnaseq_pe_yes = run_filter(
            module, tmp_path,
            "RNAseq_PE",
            RNASEQ_PE_HEADER, RNASEQ_PE_ROW,
            "yes",
            "rnaseq_pe_yes",
        )
        assert_file_absent(
            rnaseq_pe_yes,
            "RNAseq_PE + Assembly_of_ucaRNAs=yes",
            failures,
        )

        rnaseq_se_yes = run_filter(
            module, tmp_path,
            "RNAseq_SE",
            RNASEQ_SE_HEADER, RNASEQ_SE_ROW,
            "yes",
            "rnaseq_se_yes",
        )
        assert_file_absent(
            rnaseq_se_yes,
            "RNAseq_SE + Assembly_of_ucaRNAs=yes",
            failures,
        )

        ota_se_yes = run_filter(
            module, tmp_path,
            "OTA_SE",
            OTA_SE_HEADER, OTA_SE_ROW,
            "yes",
            "ota_se_yes",
        )
        assert_file_absent(
            ota_se_yes,
            "OTA_SE + Assembly_of_ucaRNAs=yes",
            failures,
        )

        ota_pe_yes = run_filter(
            module, tmp_path,
            "OTA_PE",
            OTA_PE_HEADER, OTA_PE_ROW,
            "yes",
            "ota_pe_yes",
        )
        assert_file_absent(
            ota_pe_yes,
            "OTA_PE + Assembly_of_ucaRNAs=yes",
            failures,
        )

        ata_no = run_filter(
            module, tmp_path,
            "ATA, iMARGI",
            ATA_HEADER, ATA_ROW,
            "no",
            "ata_no",
        )
        assert_file_absent(
            ata_no,
            "ATA, iMARGI + Assembly_of_ucaRNAs=no",
            failures,
        )

    if failures:
        print("FAIL: ucaRNA output scope test")
        for failure in failures:
            print("\n" + failure)
        sys.exit(1)

    print("PASS: ucaRNA output scope test")
    print("  ATA + yes:       ucaRNA ID file created")
    print("  RNAseq_PE + yes: ucaRNA ID file not created")
    print("  RNAseq_SE + yes: ucaRNA ID file not created")
    print("  OTA_SE + yes:    ucaRNA ID file not created")
    print("  OTA_PE + yes:    ucaRNA ID file not created")
    print("  ATA + no:        ucaRNA ID file not created")


if __name__ == "__main__":
    main()
