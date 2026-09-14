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
    pairtype="UU",
    rna_mapq="60",
    dna_mapq="60",
):
    return "\t".join([
        read_id,
        pairtype,

        "chr1",
        "100",
        "119",
        "+",
        "20M",
        "*",
        "*",
        "0",
        rna_mapq,

        "chr2",
        "200",
        "219",
        "+",
        "20M",
        "*",
        "*",
        "0",
        dna_mapq,
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
                f"{name:<35} "
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
            f"{name:<35} "
            f"filtered={filtered_ids!s:<30} "
            f"out={out_ids!s:<20} "
            f"{'OK' if ok else 'FAIL'}"
        )

        if not ok:
            print(
                f"    expected filtered: {expected_filtered_ids}"
            )
            print(
                f"    expected out:      {expected_out_ids}"
            )

        return ok


results = []


# ------------------------------------------------------------
# 1. Invalid r1 MAPQ between two valid records.
#    Main purpose: bad MAPQ must not stop file processing.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r1 MAPQ = *",
        [
            make_row("good_before"),
            make_row(
                "bad",
                rna_mapq="*",
            ),
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


# ------------------------------------------------------------
# 2. Invalid r2 MAPQ in a UU pair.
#    Whole pair must be rejected.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r2 MAPQ = *",
        [
            make_row(
                "bad",
                dna_mapq="*",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# ------------------------------------------------------------
# 3. Non-numeric r1 MAPQ.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r1 MAPQ = NA",
        [
            make_row(
                "bad",
                rna_mapq="NA",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# ------------------------------------------------------------
# 4. Empty r1 MAPQ.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r1 MAPQ = empty",
        [
            make_row(
                "bad",
                rna_mapq="",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# ------------------------------------------------------------
# 5. Floating-point MAPQ.
#    Current pipeline expects int(...), so this should be rejected.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r1 MAPQ = 60.0",
        [
            make_row(
                "bad",
                rna_mapq="60.0",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


# ------------------------------------------------------------
# 6. Negative MAPQ.
#    Syntactically int-compatible, but biologically/SAM-wise invalid.
# ------------------------------------------------------------

results.append(
    run_case(
        "invalid r1 MAPQ = -1",
        [
            make_row(
                "bad",
                rna_mapq="-1",
            ),
        ],
        expected_filtered_ids=[],
        expected_out_ids=[
            "bad",
        ],
    )
)


print()
print(
    f"{sum(results)}/{len(results)} tests passed"
)

if not all(results):
    sys.exit(1)