import os
import re
import sys
import tempfile
from collections import Counter

import pandas as pd
from matplotlib import pyplot as plt

from test_loader import load_functions


SCRIPT = "EditDistance_CIGAR_filter.py"

namespace = load_functions(
    SCRIPT,
    {
        "re": re,
        "Counter": Counter,
        "pd": pd,
        "plt": plt,
    }
)

namespace["save_cigar_statistics"] = lambda *args, **kwargs: None
namespace["plot_N_softClipp_or_NM"] = lambda *args, **kwargs: None

editDistance_and_CIGAR_filter = namespace["editDistance_and_CIGAR_filter"]


# ============================================================
# Headers
# ============================================================

HEADERS = {
    "ATA, not iMARGI": [
        "read_id", "ATA_pairtype",
        "rna_chr", "rna_start", "rna_end", "rna_strand",
        "rna_cigar", "rna_secondary_alignments",
        "rna_other_tags", "rna_NM", "rna_mapq",
        "dna_chr", "dna_start", "dna_end", "dna_strand",
        "dna_cigar", "dna_secondary_alignments",
        "dna_other_tags", "dna_NM", "dna_mapq",
    ],

    "ATA, iMARGI": [
        "read_id", "ATA_pairtype",
        "rna_chr", "rna_start", "rna_end", "rna_strand",
        "rna_cigar", "rna_secondary_alignments",
        "rna_other_tags", "rna_NM", "rna_mapq",
        "dna_chr", "dna_start", "dna_end", "dna_strand",
        "dna_cigar", "dna_secondary_alignments",
        "dna_other_tags", "dna_NM", "dna_mapq",
    ],

    "OTA_SE": [
        "read_id", "OTA_SE_pairtype",
        "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
        "dna1_cigar", "dna1_secondary_alignments",
        "dna1_other_tags", "dna1_NM", "dna1_mapq",
        "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
        "dna2_cigar", "dna2_secondary_alignments",
        "dna2_other_tags", "dna2_NM", "dna2_mapq",
    ],

    "OTA_PE": [
        "read_id", "OTA_PE_pairtype",
        "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
        "dna1_cigar", "dna1_secondary_alignments",
        "dna1_other_tags", "dna1_NM", "dna1_mapq",
        "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
        "dna2_cigar", "dna2_secondary_alignments",
        "dna2_other_tags", "dna2_NM", "dna2_mapq",
    ],

    "RNAseq_SE": [
        "read_id", "RNAseq_SE_pairtype",
        "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
        "rna1_cigar", "rna1_secondary_alignments",
        "rna1_other_tags", "rna1_NM", "rna1_mapq",
        "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
        "rna2_cigar", "rna2_secondary_alignments",
        "rna2_other_tags", "rna2_NM", "rna2_mapq",
    ],

    "RNAseq_PE": [
        "read_id", "RNAseq_PE_pairtype",
        "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
        "rna1_cigar", "rna1_secondary_alignments",
        "rna1_other_tags", "rna1_NM", "rna1_mapq",
        "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
        "rna2_cigar", "rna2_secondary_alignments",
        "rna2_other_tags", "rna2_NM", "rna2_mapq",
    ],
}


def make_row(
    experiment_type,
    read_id,
    pairtype,
    r1_cigar="20M",
    r2_cigar="20M",
    r1_start="100",
    r1_end="119",
    r2_start="200",
    r2_end="219",
    r1_strand="+",
    r2_strand="+",
):
    return "\t".join([
        read_id,
        pairtype,

        "chr1",
        r1_start,
        r1_end,
        r1_strand,
        r1_cigar,
        "*",
        "*",
        "0",
        "60",

        "chr1",
        r2_start,
        r2_end,
        r2_strand,
        r2_cigar,
        "*",
        "*",
        "0",
        "60",
    ])


def run_case(
    name,
    experiment_type,
    row,
    expected_filtered,
):
    with tempfile.TemporaryDirectory() as tmpdir:

        input_name = "input.tab"
        input_path = os.path.join(tmpdir, input_name)

        with open(input_path, "w") as f:
            f.write("\t".join(HEADERS[experiment_type]) + "\n")
            f.write(row + "\n")

        try:
            editDistance_and_CIGAR_filter(
                "NM + N_softClipp_bp",
                50,
                50,
                0,
                0,
                500,
                "no",
                "not explorer",
                experiment_type,
                input_name,
                tmpdir + "/",
                tmpdir + "/",
            )

        except Exception as e:
            print(
                f"{name:<48} "
                f"CRASH: {type(e).__name__}: {e}"
            )
            return False

        filtered_path = os.path.join(
            tmpdir,
            "filtered_" + input_name,
        )

        out_path = os.path.join(
            tmpdir,
            "out_" + input_name,
        )

        with open(filtered_path) as f:
            filtered_ids = [
                line.rstrip("\n").split("\t")[0]
                for line in f.readlines()[1:]
                if line.strip()
            ]

        with open(out_path) as f:
            out_ids = [
                line.rstrip("\n").split("\t")[0]
                for line in f.readlines()[1:]
                if line.strip()
            ]

        if expected_filtered:
            expected_filtered_ids = ["test"]
            expected_out_ids = []
        else:
            expected_filtered_ids = []
            expected_out_ids = ["test"]

        ok = (
            filtered_ids == expected_filtered_ids
            and out_ids == expected_out_ids
        )

        print(
            f"{name:<48} "
            f"{'OK' if ok else 'FAIL'}"
        )

        if not ok:
            print(
                f"    filtered: {filtered_ids}, "
                f"expected: {expected_filtered_ids}"
            )
            print(
                f"    out:      {out_ids}, "
                f"expected: {expected_out_ids}"
            )

        return ok


cases = []


# ============================================================
# ATA, not iMARGI
# ============================================================

cases.append((
    "ATA UU ordinary",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
    ),
    True,
))

cases.append((
    "ATA UM",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UM",
        r2_cigar="*",
    ),
    True,
))

cases.append((
    "ATA U",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "U",
        r2_cigar="*",
    ),
    True,
))

# RNA may contain one N
cases.append((
    "ATA RNA single N allowed",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="10M100N10M",
        r1_end="219",
    ),
    True,
))

# DNA may not contain N
cases.append((
    "ATA DNA single N rejected",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r2_cigar="10M100N10M",
        r2_end="319",
    ),
    False,
))

# More than one N in RNA rejected
cases.append((
    "ATA RNA multiple N rejected",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="5M50N5M50N5M",
        r1_end="214",
    ),
    False,
))


# ============================================================
# ATA, iMARGI
# ============================================================

cases.append((
    "iMARGI terminal H+S",
    "ATA, iMARGI",
    make_row(
        "ATA, iMARGI",
        "test",
        "UU",
        r1_cigar="5H10S20M",
    ),
    True,
))

cases.append((
    "iMARGI single N",
    "ATA, iMARGI",
    make_row(
        "ATA, iMARGI",
        "test",
        "UU",
        r1_cigar="10M100N10M",
        r1_end="219",
    ),
    True,
))


# ============================================================
# OTA_SE
# ============================================================

cases.append((
    "OTA_SE ordinary",
    "OTA_SE",
    make_row(
        "OTA_SE",
        "test",
        "U",
        r2_cigar="*",
    ),
    True,
))

cases.append((
    "OTA_SE N rejected",
    "OTA_SE",
    make_row(
        "OTA_SE",
        "test",
        "U",
        r1_cigar="10M100N10M",
        r1_end="219",
        r2_cigar="*",
    ),
    False,
))


# ============================================================
# OTA_PE
# ============================================================

cases.append((
    "OTA_PE UU",
    "OTA_PE",
    make_row(
        "OTA_PE",
        "test",
        "UU",
        r1_start="100",
        r1_end="119",
        r2_start="130",
        r2_end="149",
    ),
    True,
))

cases.append((
    "OTA_PE N in read1 rejected",
    "OTA_PE",
    make_row(
        "OTA_PE",
        "test",
        "UU",
        r1_cigar="10M100N10M",
        r1_start="100",
        r1_end="219",
        r2_start="230",
        r2_end="249",
    ),
    False,
))


# ============================================================
# RNAseq_SE
# ============================================================

cases.append((
    "RNAseq_SE ordinary",
    "RNAseq_SE",
    make_row(
        "RNAseq_SE",
        "test",
        "U",
        r2_cigar="*",
    ),
    True,
))

cases.append((
    "RNAseq_SE single N allowed",
    "RNAseq_SE",
    make_row(
        "RNAseq_SE",
        "test",
        "U",
        r1_cigar="10M100N10M",
        r1_end="219",
        r2_cigar="*",
    ),
    True,
))

cases.append((
    "RNAseq_SE multiple N rejected",
    "RNAseq_SE",
    make_row(
        "RNAseq_SE",
        "test",
        "U",
        r1_cigar="5M50N5M50N5M",
        r1_end="214",
        r2_cigar="*",
    ),
    False,
))


# ============================================================
# RNAseq_PE
#
# Current implementation filters/output uses rna1 only.
# ============================================================

cases.append((
    "RNAseq_PE rna1 ordinary",
    "RNAseq_PE",
    make_row(
        "RNAseq_PE",
        "test",
        "UU",
    ),
    True,
))

cases.append((
    "RNAseq_PE rna1 single N allowed",
    "RNAseq_PE",
    make_row(
        "RNAseq_PE",
        "test",
        "UU",
        r1_cigar="10M100N10M",
        r1_end="219",
    ),
    True,
))

cases.append((
    "RNAseq_PE rna1 multiple N rejected",
    "RNAseq_PE",
    make_row(
        "RNAseq_PE",
        "test",
        "UU",
        r1_cigar="5M50N5M50N5M",
        r1_end="214",
    ),
    False,
))


# ============================================================
# Other CIGAR operators
# ============================================================

cases.append((
    "Insertion accepted",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="10M2I10M",
    ),
    True,
))

cases.append((
    "Deletion accepted",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="10M2D10M",
        r1_end="121",
    ),
    True,
))

cases.append((
    "Equals and X accepted",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="10=1X9=",
    ),
    True,
))

cases.append((
    "Terminal S accepted",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="5S20M",
    ),
    True,
))

cases.append((
    "Terminal H+S accepted",
    "ATA, not iMARGI",
    make_row(
        "ATA, not iMARGI",
        "test",
        "UU",
        r1_cigar="5H10S20M",
    ),
    True,
))


# ============================================================
# Run
# ============================================================

results = []

print()
print("Experiment / pairtype / CIGAR matrix")
print("-" * 60)

for name, experiment_type, row, expected_filtered in cases:
    results.append(
        run_case(
            name,
            experiment_type,
            row,
            expected_filtered,
        )
    )

print()
print(f"{sum(results)}/{len(results)} tests passed")

if not all(results):
    sys.exit(1)