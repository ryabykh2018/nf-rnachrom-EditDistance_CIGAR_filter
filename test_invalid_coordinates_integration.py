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


HEADER = "\t".join([
    "read_id",
    "ATA_pairtype",

    "rna_chr",
    "rna_start",
    "rna_end",
    "rna_strand",
    "rna_cigar",
    "rna_secondary_alignments",
    "rna_other_tags",
    "rna_NM",
    "rna_mapq",

    "dna_chr",
    "dna_start",
    "dna_end",
    "dna_strand",
    "dna_cigar",
    "dna_secondary_alignments",
    "dna_other_tags",
    "dna_NM",
    "dna_mapq",
])


def make_row(
    read_id,
    rna_start="100",
    rna_end="119",
    dna_start="200",
    dna_end="219",
):
    return "\t".join([
        read_id,
        "UU",

        "chr1",
        rna_start,
        rna_end,
        "+",
        "20M",
        "*",
        "*",
        "0",
        "60",

        "chr2",
        dna_start,
        dna_end,
        "+",
        "20M",
        "*",
        "*",
        "0",
        "60",
    ])


def run_case(
    name,
    rows,
    expected_filtered_ids,
    expected_out_ids,
):
    with tempfile.TemporaryDirectory() as tmpdir:

        input_name = "input.tab"
        input_path = os.path.join(tmpdir, input_name)

        with open(input_path, "w") as f:
            f.write(HEADER + "\n")

            for row in rows:
                f.write(row + "\n")

        try:
            editDistance_and_CIGAR_filter(
                "NM + N_softClipp_bp",
                10,
                10,
                0,
                0,
                200,
                "no",
                "not explorer",
                "ATA, not iMARGI",
                input_name,
                tmpdir + "/",
                tmpdir + "/",
            )

        except Exception as e:
            print(
                f"{name:<40} "
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
            filtered_lines = [
                line.rstrip("\n")
                for line in f
            ]

        with open(out_path) as f:
            out_lines = [
                line.rstrip("\n")
                for line in f
            ]

        filtered_ids = [
            line.split("\t")[0]
            for line in filtered_lines[1:]
            if line
        ]

        out_ids = [
            line.split("\t")[0]
            for line in out_lines[1:]
            if line
        ]

        ok = (
            filtered_ids == expected_filtered_ids
            and out_ids == expected_out_ids
        )

        print(
            f"{name:<40} "
            f"filtered={filtered_ids!s:<30} "
            f"out={out_ids!s:<20} "
            f"{'OK' if ok else 'FAIL'}"
        )

        return ok


results = []


# 1. Bad r1 start between two good records
results.append(
    run_case(
        "invalid r1 start = *",
        [
            make_row("good_before"),
            make_row("bad", rna_start="*"),
            make_row("good_after"),
        ],
        expected_filtered_ids=[
            "good_before",
            "good_after",
        ],
        expected_out_ids=[
            "bad",
        ],
    )
)


# 2. Bad r2 end
results.append(
    run_case(
        "invalid r2 end = NA",
        [
            make_row(
                "bad",
                dna_end="NA",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# 3. Empty coordinate
results.append(
    run_case(
        "invalid r1 end = empty",
        [
            make_row(
                "bad",
                rna_end="",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# 4. start > end
results.append(
    run_case(
        "invalid r1 start > end",
        [
            make_row(
                "bad",
                rna_start="120",
                rna_end="100",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# 5. zero coordinate
results.append(
    run_case(
        "invalid r1 start = 0",
        [
            make_row(
                "bad",
                rna_start="0",
                rna_end="19",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# 6. negative coordinate
results.append(
    run_case(
        "invalid r1 start = -10",
        [
            make_row(
                "bad",
                rna_start="-10",
                rna_end="9",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


print()
print(f"{sum(results)}/{len(results)} tests passed")

if not all(results):
    sys.exit(1)