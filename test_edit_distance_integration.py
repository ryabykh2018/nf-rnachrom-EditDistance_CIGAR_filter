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


ATA_HEADER = [
    "read_id", "ATA_pairtype",
    "rna_chr", "rna_start", "rna_end", "rna_strand",
    "rna_cigar", "rna_NM", "rna_mapq",
    "dna_chr", "dna_start", "dna_end", "dna_strand",
    "dna_cigar", "dna_NM", "dna_mapq",
    "rna_secondary_alignments", "dna_secondary_alignments",
    "rna_other_tags", "dna_other_tags"
]


def write_tsv(path, rows):
    with open(path, "w") as f:
        f.write("\t".join(ATA_HEADER) + "\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


def read_ids(path):
    with open(path) as f:
        next(f)
        return {
            line.rstrip("\n").split("\t")[0]
            for line in f
            if line.strip()
        }


def run_case(name, experiment_type, rows, expected_filtered, expected_out):
    with tempfile.TemporaryDirectory() as tmp:
        input_dir = os.path.join(tmp, "input")
        output_dir = os.path.join(tmp, "output")

        os.makedirs(input_dir)
        os.makedirs(output_dir)

        filename = name + ".tsv"

        write_tsv(
            os.path.join(input_dir, filename),
            rows
        )

        editDistance_and_CIGAR_filter(
            "NM + N_softClipp_bp",
            3,      # r1 threshold
            3,      # r2 threshold
            0,
            0,
            200,
            "no",
            "not explorer",
            experiment_type,
            filename,
            input_dir + "/",
            output_dir + "/"
        )

        filtered = read_ids(
            os.path.join(output_dir, "filtered_" + filename)
        )

        out = read_ids(
            os.path.join(output_dir, "out_" + filename)
        )

        ok = (
            filtered == set(expected_filtered)
            and out == set(expected_out)
        )

        print(f"\n{name}")
        print("-" * len(name))
        print("filtered expected:", sorted(expected_filtered))
        print("filtered got:     ", sorted(filtered))
        print("out expected:     ", sorted(expected_out))
        print("out got:          ", sorted(out))
        print("RESULT:", "OK" if ok else "FAIL")

        return ok


def row(read_id, strand, cigar, nm):
    return [
        read_id, "UU",

        # RNA
        "chr1", 100, 119, strand,
        cigar, nm, 60,

        # DNA — deliberately boring and always passing
        "chr1", 300, 319, "+",
        "20M", 0, 60,

        "*", "*", "*", "*"
    ]


# ------------------------------------------------------------
# ATA, not iMARGI
#
# threshold = 3
# all clipping counts
# ------------------------------------------------------------

ATA_NORMAL_ROWS = [
    row("normal_20M", "+", "20M", 3),             # 3 -> PASS
    row("normal_2S20M", "+", "2S20M", 1),         # 1+2=3 -> PASS
    row("normal_3S20M", "+", "3S20M", 1),         # 1+3=4 -> FAIL
    row("normal_20M2S", "+", "20M2S", 1),         # 1+2=3 -> PASS
    row("normal_double", "+", "1S20M2S", 1),       # 1+1+2=4 -> FAIL
    row("normal_HS", "+", "1H2S20M", 0),           # 3 -> PASS
    row("normal_HS_fail", "+", "2H2S20M", 0),      # 4 -> FAIL
]


# ------------------------------------------------------------
# ATA, iMARGI
#
# + strand -> RIGHT terminal clipping ignored
# - strand -> LEFT terminal clipping ignored
# ------------------------------------------------------------

IMARGI_ROWS = [
    # +
    row("imargi_plus_left_pass", "+", "2S20M", 1),        # 1+2=3 PASS
    row("imargi_plus_left_fail", "+", "3S20M", 1),        # 1+3=4 FAIL

    row("imargi_plus_right_ignored", "+", "20M10S", 3),    # 3 PASS
    row("imargi_plus_double", "+", "2S20M10S", 1),        # right ignored -> 1+2=3 PASS

    # -
    row("imargi_minus_left_ignored", "-", "10S20M", 3),    # 3 PASS
    row("imargi_minus_right_pass", "-", "20M2S", 1),       # 1+2=3 PASS
    row("imargi_minus_right_fail", "-", "20M3S", 1),       # 1+3=4 FAIL

    # H + S groups
    row("imargi_plus_HS_ignored", "+", "20M2S2H", 3),      # right group ignored
    row("imargi_minus_HS_counted", "-", "20M2S1H", 0),     # right group = 3 PASS
    row("imargi_minus_HS_fail", "-", "20M2S2H", 0),        # right group = 4 FAIL
]


results = []

results.append(
    run_case(
        "ATA_normal_ED",
        "ATA, not iMARGI",
        ATA_NORMAL_ROWS,
        expected_filtered={
            "normal_20M",
            "normal_2S20M",
            "normal_20M2S",
            "normal_HS",
        },
        expected_out={
            "normal_3S20M",
            "normal_double",
            "normal_HS_fail",
        }
    )
)

results.append(
    run_case(
        "ATA_iMARGI_ED",
        "ATA, iMARGI",
        IMARGI_ROWS,
        expected_filtered={
            "imargi_plus_left_pass",
            "imargi_plus_right_ignored",
            "imargi_plus_double",
            "imargi_minus_left_ignored",
            "imargi_minus_right_pass",
            "imargi_plus_HS_ignored",
            "imargi_minus_HS_counted",
        },
        expected_out={
            "imargi_plus_left_fail",
            "imargi_minus_right_fail",
            "imargi_minus_HS_fail",
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