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


HEADER = [
    "read_id", "ATA_pairtype",
    "rna_chr", "rna_start", "rna_end", "rna_strand",
    "rna_cigar", "rna_NM", "rna_mapq",
    "dna_chr", "dna_start", "dna_end", "dna_strand",
    "dna_cigar", "dna_NM", "dna_mapq",
    "rna_secondary_alignments", "dna_secondary_alignments",
    "rna_other_tags", "dna_other_tags"
]


def make_row(read_id, cigar, start, end, strand="+"):
    return [
        read_id, "UU",

        # RNA
        "chr1", start, end, strand,
        cigar, 0, 60,

        # DNA: simple passing mate
        "chr1", 1000, 1019, "+",
        "20M", 0, 60,

        "*", "*", "*", "*"
    ]


tests = [
    # read_id, cigar, input start, input end, expected start, expected end

    (
        "left_longer",
        "30M100N28M",
        100,
        257,
        100,
        129
    ),

    (
        "right_longer",
        "3M100N4M",
        100,
        206,
        203,
        206
    ),

    (
        "left_clip_left_longer",
        "2S30M100N28M",
        100,
        257,
        100,
        129
    ),

    (
        "right_clip_right_longer",
        "3M100N4M2S",
        100,
        206,
        203,
        206
    ),

    # I does not consume reference
    (
        "left_with_I",
        "10M3I10M100N19M",
        100,
        238,
        100,
        119
    ),

    # D consumes reference
    (
        "left_with_D",
        "10M3D10M100N19M",
        100,
        241,
        100,
        122
    ),

    # = and X behave as aligned/reference-consuming
    (
        "left_with_equal_X",
        "10=1X10=100N19M",
        100,
        239,
        100,
        120
    ),
]


with tempfile.TemporaryDirectory() as tmp:
    input_dir = os.path.join(tmp, "input")
    output_dir = os.path.join(tmp, "output")

    os.makedirs(input_dir)
    os.makedirs(output_dir)

    filename = "single_N_coords.tsv"

    rows = [
        make_row(
            read_id,
            cigar,
            start,
            end
        )
        for read_id, cigar, start, end, exp_start, exp_end in tests
    ]

    input_path = os.path.join(input_dir, filename)

    with open(input_path, "w") as f:
        f.write("\t".join(HEADER) + "\n")

        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


    editDistance_and_CIGAR_filter(
        "NM + N_softClipp_bp",
        100,    # r1 ED threshold
        100,    # r2 ED threshold
        0,
        0,
        200,
        "no",
        "not explorer",
        "ATA, not iMARGI",
        filename,
        input_dir + "/",
        output_dir + "/"
    )


    filtered_path = os.path.join(
        output_dir,
        "filtered_" + filename
    )

    with open(filtered_path) as f:
        header = f.readline().rstrip("\n").split("\t")
        header_dict = {
            name: i
            for i, name in enumerate(header)
        }

        results = {}

        for line in f:
            fields = line.rstrip("\n").split("\t")

            read_id = fields[header_dict["read_id"]]

            results[read_id] = (
                int(fields[header_dict["rna_start"]]),
                int(fields[header_dict["rna_end"]])
            )


print(
    f"{'read_id':<25} "
    f"{'expected':<15} "
    f"{'got':<15} "
    f"RESULT"
)

print("-" * 70)

n_failed = 0

for (
    read_id,
    cigar,
    start,
    end,
    exp_start,
    exp_end
) in tests:

    expected = (exp_start, exp_end)
    got = results.get(read_id)

    ok = got == expected

    if not ok:
        n_failed += 1

    print(
        f"{read_id:<25} "
        f"{str(expected):<15} "
        f"{str(got):<15} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if not all(results):
    sys.exit(1)