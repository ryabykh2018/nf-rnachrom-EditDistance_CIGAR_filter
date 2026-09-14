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

parse_cigar_tokens = namespace["parse_cigar_tokens"]


tests = [
    # CIGAR, expected valid

    ("20M", True),
    ("2S20M", True),
    ("1H2S20M", True),
    ("20M2S1H", True),
    ("1H2S20M2S1H", True),

    ("10M2I8M", True),
    ("10M2D10M", True),
    ("10=2X8=", True),
    ("10M100N10M", True),

    ("2S10M2I5M3D4=1X100N8M1S", True),

    # invalid / unsupported
    ("*", False),
    ("", False),

    ("10MFOO5S", False),
    ("10M5", False),
    ("M10", False),
    ("10M2Q5M", False),

    # P пока intentionally unsupported
    ("10M2P10M", False),

    # zero-length operation
    ("0M", False),
    ("10M0S", False),
]


print(
    f"{'CIGAR':<30} "
    f"{'expected valid':<15} "
    f"{'got valid':<12} "
    f"RESULT"
)

print("-" * 75)

n_failed = 0

for cigar, expected_valid in tests:

    result = parse_cigar_tokens(cigar)
    got_valid = result is not None

    ok = got_valid == expected_valid

    if not ok:
        n_failed += 1

    print(
        f"{repr(cigar):<30} "
        f"{str(expected_valid):<15} "
        f"{str(got_valid):<12} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)