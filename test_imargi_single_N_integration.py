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


def make_row(read_id, cigar, start, end, strand, nm):
    return [
        read_id, "UU",

        "chr1", start, end, strand,
        cigar, nm, 60,

        "chr1", 1000, 1019, "+",
        "20M", 0, 60,

        "*", "*", "*", "*"
    ]


tests = [
    # read_id, strand, cigar, start, end, NM,
    # expected start, expected end, expected pass

    # + strand:
    # right clip = 3' technical end -> ignored for ED

    (
        "plus_left_longer_right_clip",
        "+",
        "30M100N28M5S",
        100, 257,
        3,
        100, 129,
        True
    ),

    (
        "plus_left_longer_left_clip",
        "+",
        "2S30M100N28M",
        100, 257,
        1,
        100, 129,
        True
    ),

    # left clip counts: NM 2 + 2S = 4 -> fail at threshold 3
    (
        "plus_left_clip_fail",
        "+",
        "2S30M100N28M",
        100, 257,
        2,
        100, 129,
        False
    ),

    # right block longer; right clip still ignored for ED
    (
        "plus_right_longer_right_clip",
        "+",
        "3M100N4M5S",
        100, 206,
        3,
        203, 206,
        True
    ),

    # - strand:
    # left clip = 3' technical end -> ignored for ED

    (
        "minus_left_longer_left_clip",
        "-",
        "5S30M100N28M",
        100, 257,
        3,
        100, 129,
        True
    ),

    (
        "minus_right_longer_right_clip",
        "-",
        "3M100N4M2S",
        100, 206,
        1,
        203, 206,
        True
    ),

    # right clip counts: NM 2 + 2S = 4 -> fail
    (
        "minus_right_clip_fail",
        "-",
        "3M100N4M2S",
        100, 206,
        2,
        203, 206,
        False
    ),

    # H+S technical clipping group
    (
        "plus_HS_technical",
        "+",
        "30M100N28M2S1H",
        100, 257,
        3,
        100, 129,
        True
    ),

    (
        "minus_HS_technical",
        "-",
        "1H2S30M100N28M",
        100, 257,
        3,
        100, 129,
        True
    ),
]


with tempfile.TemporaryDirectory() as tmp:

    input_dir = os.path.join(tmp, "input")
    output_dir = os.path.join(tmp, "output")

    os.makedirs(input_dir)
    os.makedirs(output_dir)

    filename = "imargi_single_N.tsv"

    rows = [
        make_row(
            read_id,
            cigar,
            start,
            end,
            strand,
            nm
        )
        for (
            read_id,
            strand,
            cigar,
            start,
            end,
            nm,
            exp_start,
            exp_end,
            exp_pass
        ) in tests
    ]

    with open(os.path.join(input_dir, filename), "w") as f:
        f.write("\t".join(HEADER) + "\n")

        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")


    editDistance_and_CIGAR_filter(
        "NM + N_softClipp_bp",
        3,      # r1 ED threshold
        3,      # r2 ED threshold
        0,
        0,
        200,
        "no",
        "not explorer",
        "ATA, iMARGI",
        filename,
        input_dir + "/",
        output_dir + "/"
    )


    def read_table(path):
        result = {}

        with open(path) as f:
            header = f.readline().rstrip("\n").split("\t")
            hd = {name: i for i, name in enumerate(header)}

            for line in f:
                fields = line.rstrip("\n").split("\t")

                read_id = fields[hd["read_id"]]

                result[read_id] = (
                    int(fields[hd["rna_start"]]),
                    int(fields[hd["rna_end"]])
                )

        return result


    filtered = read_table(
        os.path.join(output_dir, "filtered_" + filename)
    )

    out = read_table(
        os.path.join(output_dir, "out_" + filename)
    )


print(
    f"{'read_id':<32} "
    f"{'pass exp':<9} "
    f"{'coords exp':<15} "
    f"{'where got':<10} "
    f"{'coords got':<15} "
    f"RESULT"
)

print("-" * 100)

n_failed = 0

for (
    read_id,
    strand,
    cigar,
    start,
    end,
    nm,
    exp_start,
    exp_end,
    exp_pass
) in tests:

    expected_coords = (exp_start, exp_end)

    if read_id in filtered:
        where = "filtered"
        got_coords = filtered[read_id]
        got_pass = True

    elif read_id in out:
        where = "out"
        got_coords = out[read_id]
        got_pass = False

    else:
        where = "missing"
        got_coords = None
        got_pass = None

    # For rejected rows output keeps original coordinates,
    # so only compare corrected coordinates for passed rows.
    coords_ok = (
        got_coords == expected_coords
        if exp_pass
        else True
    )

    ok = (
        got_pass == exp_pass
        and coords_ok
    )

    if not ok:
        n_failed += 1

    print(
        f"{read_id:<32} "
        f"{str(exp_pass):<9} "
        f"{str(expected_coords):<15} "
        f"{where:<10} "
        f"{str(got_coords):<15} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)