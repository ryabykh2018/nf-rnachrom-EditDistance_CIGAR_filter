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
calculate_edit_distance = namespace["calculate_edit_distance"]


# ------------------------------------------------------------
# Tests
#
# experiment, strand, CIGAR, NM,
# expected clipping, expected final ED
# ------------------------------------------------------------

tests = [

    # ========================================================
    # Ordinary mode: all terminal S/H count
    # ========================================================

    ("ATA, not iMARGI", "+", "20M",        2, 0, 2),
    ("ATA, not iMARGI", "+", "2S20M",      2, 2, 4),
    ("ATA, not iMARGI", "+", "20M2S",      2, 2, 4),
    ("ATA, not iMARGI", "+", "2S20M3S",    2, 5, 7),

    ("ATA, not iMARGI", "-", "2S20M3S",    2, 5, 7),

    # H + S on same side
    ("ATA, not iMARGI", "+", "1H2S20M",    1, 3, 4),
    ("ATA, not iMARGI", "+", "20M2S1H",    1, 3, 4),
    ("ATA, not iMARGI", "+", "1H2S20M2S1H", 1, 6, 7),

    # ========================================================
    # ATA iMARGI
    #
    # + strand -> ignore RIGHT terminal clipping
    # - strand -> ignore LEFT terminal clipping
    # ========================================================

    ("ATA, iMARGI", "+", "2S20M",          2, 2, 4),
    ("ATA, iMARGI", "+", "20M2S",          2, 0, 2),
    ("ATA, iMARGI", "+", "2S20M3S",        2, 2, 4),

    ("ATA, iMARGI", "-", "2S20M",          2, 0, 2),
    ("ATA, iMARGI", "-", "20M2S",          2, 2, 4),
    ("ATA, iMARGI", "-", "2S20M3S",        2, 3, 5),

    # H + S terminal group
    ("ATA, iMARGI", "+", "1H2S20M",        1, 3, 4),
    ("ATA, iMARGI", "+", "20M2S1H",        1, 0, 1),

    ("ATA, iMARGI", "-", "1H2S20M",        1, 0, 1),
    ("ATA, iMARGI", "-", "20M2S1H",        1, 3, 4),

    ("ATA, iMARGI", "+", "1H2S20M2S1H",    1, 3, 4),
    ("ATA, iMARGI", "-", "1H2S20M2S1H",    1, 3, 4),

    # ========================================================
    # I / D / = / X
    #
    # These must not add anything beyond NM.
    # ========================================================

    ("ATA, not iMARGI", "+", "10M2I8M",       2, 0, 2),
    ("ATA, not iMARGI", "+", "10M2D10M",      2, 0, 2),
    ("ATA, not iMARGI", "+", "10=2X8=",       2, 0, 2),

    ("ATA, not iMARGI", "+", "2S10M2I8M",     2, 2, 4),
    ("ATA, not iMARGI", "+", "10M2D10M3S",    2, 3, 5),
    ("ATA, not iMARGI", "+", "1H2S10=2X8=",   2, 3, 5),

    # ========================================================
    # OTA / RNAseq:
    # clipping rules should be ordinary, not iMARGI-specific
    # ========================================================

    ("OTA_SE", "+", "2S20M3S",       1, 5, 6),
    ("OTA_PE", "-", "1H2S20M",       1, 3, 4),

    ("RNAseq_SE", "+", "2S20M3S",    1, 5, 6),
    ("RNAseq_PE", "-", "20M2S1H",    1, 3, 4),
]


print(
    f"{'experiment':<18} "
    f"{'str':<4} "
    f"{'CIGAR':<18} "
    f"{'clip exp':<9} "
    f"{'clip got':<9} "
    f"{'ED exp':<7} "
    f"{'ED got':<7} "
    f"RESULT"
)

print("-" * 90)

n_failed = 0

for experiment, strand, cigar, nm, expected_clip, expected_ed in tests:

    cigar_type, clipping, clipping_type = CIGAR_field_classifier(
        cigar,
        experiment,
        strand
    )

    got_clip = sum(clipping)

    got_ed = calculate_edit_distance(
        nm,
        clipping,
        "NM + N_softClipp_bp"
    )

    ok = (
        got_clip == expected_clip
        and got_ed == expected_ed
    )

    if not ok:
        n_failed += 1

    print(
        f"{experiment:<18} "
        f"{strand:<4} "
        f"{cigar:<18} "
        f"{expected_clip:<9} "
        f"{got_clip:<9} "
        f"{expected_ed:<7} "
        f"{got_ed:<7} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)