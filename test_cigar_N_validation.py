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

validate_cigar_N = namespace["validate_cigar_N"]


tests = [
    # experiment_type, read_role, CIGAR, expected

    # ATA, not iMARGI
    ("ATA, not iMARGI", "rna", "20M", True),
    ("ATA, not iMARGI", "rna", "10M100N10M", True),
    ("ATA, not iMARGI", "rna", "5M100N5M200N10M", False),

    ("ATA, not iMARGI", "dna", "20M", True),
    ("ATA, not iMARGI", "dna", "10M100N10M", False),
    ("ATA, not iMARGI", "dna", "5M100N5M200N10M", False),

    # ATA, iMARGI
    ("ATA, iMARGI", "rna", "20M", True),
    ("ATA, iMARGI", "rna", "10M100N10M", True),
    ("ATA, iMARGI", "rna", "5M100N5M200N10M", False),

    ("ATA, iMARGI", "dna", "20M", True),
    ("ATA, iMARGI", "dna", "10M100N10M", False),
    ("ATA, iMARGI", "dna", "5M100N5M200N10M", False),

    # OTA_SE
    ("OTA_SE", "dna", "20M", True),
    ("OTA_SE", "dna", "10M100N10M", False),
    ("OTA_SE", "dna", "5M100N5M200N10M", False),

    # OTA_PE
    ("OTA_PE", "dna1", "20M", True),
    ("OTA_PE", "dna1", "10M100N10M", False),

    ("OTA_PE", "dna2", "20M", True),
    ("OTA_PE", "dna2", "10M100N10M", False),

    # RNAseq_SE
    ("RNAseq_SE", "rna", "20M", True),
    ("RNAseq_SE", "rna", "10M100N10M", True),
    ("RNAseq_SE", "rna", "5M100N5M200N10M", False),

    # RNAseq_PE
    ("RNAseq_PE", "rna1", "20M", True),
    ("RNAseq_PE", "rna1", "10M100N10M", True),
    ("RNAseq_PE", "rna1", "5M100N5M200N10M", False),
]


print(
    f"{'experiment':<18} "
    f"{'role':<6} "
    f"{'CIGAR':<22} "
    f"{'expected':<9} "
    f"{'got':<6} "
    f"RESULT"
)

print("-" * 80)

n_failed = 0

for experiment_type, read_role, cigar, expected in tests:
    got = validate_cigar_N(cigar, experiment_type, read_role)

    ok = got == expected

    if not ok:
        n_failed += 1

    print(
        f"{experiment_type:<18} "
        f"{read_role:<6} "
        f"{cigar:<22} "
        f"{str(expected):<9} "
        f"{str(got):<6} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)