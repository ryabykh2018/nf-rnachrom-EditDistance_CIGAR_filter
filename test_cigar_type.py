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


tests = [
    # experiment, strand, cigar, expected_cigar_type

    ("ATA, not iMARGI", "+", "20M", "M1"),
    ("ATA, not iMARGI", "+", "2S20M", "M1S1"),
    ("ATA, not iMARGI", "+", "2S20M3S", "M1S2"),

    ("ATA, not iMARGI", "+", "1H2S20M", "H1M1S1"),
    ("ATA, not iMARGI", "+", "20M2S1H", "H1M1S1"),
    ("ATA, not iMARGI", "+", "1H2S20M2S1H", "H2M1S2"),

    ("ATA, not iMARGI", "+", "10M2I8M", "I1M2"),
    ("ATA, not iMARGI", "+", "10M2D10M", "D1M2"),
    ("ATA, not iMARGI", "+", "10=2X8=", "=2X1"),

    ("ATA, not iMARGI", "+", "2S10M2I8M3S", "I1M2S2"),
    ("ATA, not iMARGI", "+", "10M2D5M100N5M", "D1M3N1"),

    # iMARGI: technical-end clipping is removed before cigar_type
    ("ATA, iMARGI", "+", "2S20M3S", "M1S1"),
    ("ATA, iMARGI", "-", "2S20M3S", "M1S1"),

    ("ATA, iMARGI", "+", "1H2S20M2S1H", "H1M1S1"),
    ("ATA, iMARGI", "-", "1H2S20M2S1H", "H1M1S1"),
]


print(
    f"{'experiment':<18} "
    f"{'str':<4} "
    f"{'CIGAR':<22} "
    f"{'expected':<12} "
    f"{'got':<12} "
    f"RESULT"
)

print("-" * 80)

n_failed = 0

for experiment, strand, cigar, expected in tests:

    cigar_type, clipping, clipping_type = CIGAR_field_classifier(
        cigar,
        experiment,
        strand
    )

    ok = cigar_type == expected

    if not ok:
        n_failed += 1

    print(
        f"{experiment:<18} "
        f"{strand:<4} "
        f"{cigar:<22} "
        f"{expected:<12} "
        f"{cigar_type:<12} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)