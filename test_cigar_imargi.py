import re
import sys
from collections import Counter

from test_loader import load_functions


SCRIPT = "EditDistance_CIGAR_filter.py"

namespace = load_functions(
    SCRIPT,
    {
        "re": re,
        "Counter": Counter,
    }
)

CIGAR_field_classifier = namespace["CIGAR_field_classifier"]
coordinate_trimmer_by_cigar = namespace["coordinate_trimmer_by_cigar"]


# ------------------------------------------------------------
# Тесты
# ------------------------------------------------------------

START = 100
END = 119
EDIT_DIST_TYPE = "NM + N_softClipp_bp"
EXPERIMENT = "ATA, iMARGI"


tests = [
    # strand, CIGAR, expected clipping penalty

    ("+", "20M",       0),
    ("+", "2S20M",     2),
    ("+", "20M2S",     0),
    ("+", "2S20M2S",   2),

    ("-", "20M",       0),
    ("-", "2S20M",     0),
    ("-", "20M2S",     2),
    ("-", "2S20M2S",   2),

    # H
    ("+", "2H20M",     2),
    ("+", "20M2H",     0),

    ("-", "2H20M",     0),
    ("-", "20M2H",     2),

    # H + S на одном конце
    ("+", "1H2S20M",   3),
    ("+", "20M2S1H",   0),

    ("-", "1H2S20M",   0),
    ("-", "20M2S1H",   3),
]


print(
    f"{'strand':<7} "
    f"{'CIGAR':<12} "
    f"{'clip exp':<10} "
    f"{'clip got':<10} "
    f"{'coord exp':<12} "
    f"{'coord got':<12} "
    f"RESULT"
)

print("-" * 80)

n_failed = 0

for strand, cigar, expected_clip in tests:

    cigar_type, clipping, clipping_type = CIGAR_field_classifier(
        cigar,
        EXPERIMENT,
        strand
    )

    actual_clip = sum(clipping)

    start_new, end_new = coordinate_trimmer_by_cigar(
        EDIT_DIST_TYPE,
        cigar,
        START,
        END,
        strand,
        EXPERIMENT
    )

    expected_coords = (START, END)
    actual_coords = (start_new, end_new)

    clip_ok = actual_clip == expected_clip
    coords_ok = actual_coords == expected_coords

    ok = clip_ok and coords_ok

    if not ok:
        n_failed += 1

    print(
        f"{strand:<7} "
        f"{cigar:<12} "
        f"{expected_clip:<10} "
        f"{actual_clip:<10} "
        f"{str(expected_coords):<12} "
        f"{str(actual_coords):<12} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)