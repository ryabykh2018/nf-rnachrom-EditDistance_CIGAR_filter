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

# Statistics and plots are irrelevant for integration tests
namespace["save_cigar_statistics"] = lambda *args, **kwargs: None
namespace["plot_N_softClipp_or_NM"] = lambda *args, **kwargs: None

editDistance_and_CIGAR_filter = namespace["editDistance_and_CIGAR_filter"]



# ============================================================
# Helpers
# ============================================================

def write_tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")

        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def read_ids(path):
    with open(path) as f:
        next(f)  # header
        return {
            line.rstrip("\n").split("\t")[0]
            for line in f
            if line.strip()
        }

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

def run_case(
    name,
    experiment_type,
    header,
    rows,
    expected_filtered,
    expected_out,
    coordinate_check=None
):
    with tempfile.TemporaryDirectory() as tmp:

        input_dir = os.path.join(tmp, "input")
        output_dir = os.path.join(tmp, "output")

        os.makedirs(input_dir)
        os.makedirs(output_dir)

        filename = f"{name}.tsv"

        write_tsv(
            os.path.join(input_dir, filename),
            header,
            rows
        )

        editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",  # edit_dist_type
            100,                    # r1 edit distance threshold
            100,                    # r2 edit distance threshold
            0,                      # r1 MAPQ threshold
            0,                      # r2 MAPQ threshold
            200,                    # OTA_PE distance threshold
            "no",                   # Assembly_of_ucaRNAs
            "not explorer",         # mode
            experiment_type,
            filename,
            input_dir + "/",
            output_dir + "/"
        )

        filtered_path = os.path.join(
            output_dir,
            "filtered_" + filename
        )

        out_path = os.path.join(
            output_dir,
            "out_" + filename
        )

        filtered = read_ids(filtered_path)
        out = read_ids(out_path)

        ok = (
            filtered == set(expected_filtered)
            and out == set(expected_out)
        )

        # Optionally verify that a rejected contact preserved
        # its original genomic coordinates.
        if coordinate_check is not None:
            row = read_row_by_id(
                out_path,
                coordinate_check["read_id"]
            )

            if row is None:
                ok = False
                print(
                    "Coordinate check: FAIL -",
                    coordinate_check["read_id"],
                    "was not found in out_*"
                )
            else:
                start_ok = (
                    row[coordinate_check["start_column"]]
                    == str(coordinate_check["expected_start"])
                )

                end_ok = (
                    row[coordinate_check["end_column"]]
                    == str(coordinate_check["expected_end"])
                )

                if not (start_ok and end_ok):
                    ok = False

                print(
                    "Coordinate check:",
                    f'{row[coordinate_check["start_column"]]}-'
                    f'{row[coordinate_check["end_column"]]}',
                    "expected:",
                    f'{coordinate_check["expected_start"]}-'
                    f'{coordinate_check["expected_end"]}'
                )

        print(f"\n{name}")
        print("-" * len(name))
        print("filtered expected:", sorted(expected_filtered))
        print("filtered got:     ", sorted(filtered))
        print("out expected:     ", sorted(expected_out))
        print("out got:          ", sorted(out))
        print("RESULT:", "OK" if ok else "FAIL")

        return ok


# ============================================================
# Headers
# ============================================================

ATA_HEADER = [
    "read_id", "ATA_pairtype",
    "rna_chr", "rna_start", "rna_end", "rna_strand",
    "rna_cigar", "rna_NM", "rna_mapq",
    "dna_chr", "dna_start", "dna_end", "dna_strand",
    "dna_cigar", "dna_NM", "dna_mapq",
    "rna_secondary_alignments", "dna_secondary_alignments",
    "rna_other_tags", "dna_other_tags"
]


OTA_SE_HEADER = [
    "read_id", "OTA_SE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags"
]


OTA_PE_HEADER = [
    "read_id", "OTA_PE_pairtype",
    "dna1_chr", "dna1_start", "dna1_end", "dna1_strand",
    "dna1_cigar", "dna1_NM", "dna1_mapq",
    "dna2_chr", "dna2_start", "dna2_end", "dna2_strand",
    "dna2_cigar", "dna2_NM", "dna2_mapq",
    "dna1_secondary_alignments", "dna2_secondary_alignments",
    "dna1_other_tags", "dna2_other_tags"
]


RNASEQ_SE_HEADER = [
    "read_id", "RNAseq_SE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags"
]


RNASEQ_PE_HEADER = [
    "read_id", "RNAseq_PE_pairtype",
    "rna1_chr", "rna1_start", "rna1_end", "rna1_strand",
    "rna1_cigar", "rna1_NM", "rna1_mapq",
    "rna2_chr", "rna2_start", "rna2_end", "rna2_strand",
    "rna2_cigar", "rna2_NM", "rna2_mapq",
    "rna1_secondary_alignments", "rna2_secondary_alignments",
    "rna1_other_tags", "rna2_other_tags"
]


# ============================================================
# 1. ATA
#
# RNA:
#   0 N -> pass
#   1 N -> pass
#   2 N -> reject
#
# DNA:
#   1 N -> reject
#
# UM:
#   only RNA is checked
# ============================================================

ATA_ROWS = [

    # no N anywhere -> PASS
    [
        "ATA_noN", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 200, 219, "+", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    # RNA has one N -> PASS
    [
        "ATA_RNA_1N", "UU",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "chr1", 300, 319, "+", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    # RNA has two N -> REJECT
    [
        "ATA_RNA_2N", "UU",
        "chr1", 100, 419, "+", "5M100N5M200N10M", 0, 60,
        "chr1", 500, 519, "+", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    # DNA has one N -> REJECT
    [
        "ATA_DNA_1N", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 300, 419, "+", "10M100N10M", 0, 60,
        "*", "*", "*", "*"
    ],

    # UM: RNA has one N -> PASS.
    # DNA side is multimapper and is not filtered as r2.
    [
        "ATA_UM_RNA_1N", "UM",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],

    # Rejected contacts must preserve the original genomic coordinates in out_*
    [
        "test_out_original_coords", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 200, 318, "+", "9M100N10M", 0, 60,
        "*", "*", "*", "*"
    ],
]


# ============================================================
# 2. OTA_SE
#
# Any N -> reject
# ============================================================

OTA_SE_ROWS = [
    [
        "OTA_SE_noN", "U",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],

    [
        "OTA_SE_1N", "U",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],
]


# ============================================================
# 3. OTA_PE
#
# Both reads must be N-free.
# One bad read -> whole pair rejected.
#
# Keep mates close enough to satisfy distance threshold.
# ============================================================

OTA_PE_ROWS = [
    [
        "OTA_PE_noN", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 130, 149, "-", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    [
        "OTA_PE_dna1_1N", "UU",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "chr1", 220, 239, "-", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    [
        "OTA_PE_dna2_1N", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 120, 239, "-", "10M100N10M", 0, 60,
        "*", "*", "*", "*"
    ],
]


# ============================================================
# 4. RNAseq_SE
#
# 0/1 N -> pass
# 2 N -> reject
# ============================================================

RNASEQ_SE_ROWS = [
    [
        "RNAseq_SE_noN", "U",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],

    [
        "RNAseq_SE_1N", "U",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],

    [
        "RNAseq_SE_2N", "U",
        "chr1", 100, 419, "+", "5M100N5M200N10M", 0, 60,
        "*", "*", "*", "*", "*", "*", "*",
        "*", "*", "*", "*"
    ],
]


# ============================================================
# 5. RNAseq_PE
#
# Current filter works only with rna1.
# Therefore N-rules should be applied to rna1.
# ============================================================

RNASEQ_PE_ROWS = [
    [
        "RNAseq_PE_r1_noN", "UU",
        "chr1", 100, 119, "+", "20M", 0, 60,
        "chr1", 300, 319, "-", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    [
        "RNAseq_PE_r1_1N", "UU",
        "chr1", 100, 219, "+", "10M100N10M", 0, 60,
        "chr1", 300, 319, "-", "20M", 0, 60,
        "*", "*", "*", "*"
    ],

    [
        "RNAseq_PE_r1_2N", "UU",
        "chr1", 100, 419, "+", "5M100N5M200N10M", 0, 60,
        "chr1", 500, 519, "-", "20M", 0, 60,
        "*", "*", "*", "*"
    ],
]


# ============================================================
# Run
# ============================================================

results = []

results.append(
    run_case(
        "ATA_N_validation",
        "ATA, not iMARGI",
        ATA_HEADER,
        ATA_ROWS,
        expected_filtered={
            "ATA_noN",
            "ATA_RNA_1N",
            "ATA_UM_RNA_1N",
        },
        expected_out={
            "ATA_RNA_2N",
            "ATA_DNA_1N",
            "test_out_original_coords",
        },
        coordinate_check={
            "read_id": "test_out_original_coords",
            "start_column": "dna_start",
            "end_column": "dna_end",
            "expected_start": 200,
            "expected_end": 318,
        }
    )
)

results.append(
    run_case(
        "OTA_SE_N_validation",
        "OTA_SE",
        OTA_SE_HEADER,
        OTA_SE_ROWS,
        expected_filtered={
            "OTA_SE_noN",
        },
        expected_out={
            "OTA_SE_1N",
        }
    )
)

results.append(
    run_case(
        "OTA_PE_N_validation",
        "OTA_PE",
        OTA_PE_HEADER,
        OTA_PE_ROWS,
        expected_filtered={
            "OTA_PE_noN",
        },
        expected_out={
            "OTA_PE_dna1_1N",
            "OTA_PE_dna2_1N",
        }
    )
)

results.append(
    run_case(
        "RNAseq_SE_N_validation",
        "RNAseq_SE",
        RNASEQ_SE_HEADER,
        RNASEQ_SE_ROWS,
        expected_filtered={
            "RNAseq_SE_noN",
            "RNAseq_SE_1N",
        },
        expected_out={
            "RNAseq_SE_2N",
        }
    )
)

results.append(
    run_case(
        "RNAseq_PE_N_validation",
        "RNAseq_PE",
        RNASEQ_PE_HEADER,
        RNASEQ_PE_ROWS,
        expected_filtered={
            "RNAseq_PE_r1_noN",
            "RNAseq_PE_r1_1N",
        },
        expected_out={
            "RNAseq_PE_r1_2N",
        }
    )
)


print("\n" + "=" * 60)

n_ok = sum(results)

print(f"Integration tests: {n_ok}/{len(results)} OK")

if n_ok == len(results):
    print("ALL OK")
else:
    print("SOME TESTS FAILED")

if not all(results):
    sys.exit(1)