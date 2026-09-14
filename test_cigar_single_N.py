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

coordinate_trimmer_by_cigar = namespace["coordinate_trimmer_by_cigar"]


EDIT_DIST_TYPE = "NM + N_softClipp_bp"
EXPERIMENT = "ATA, not iMARGI"


tests = [
    # cigar, start, end, strand, expected_start, expected_end

    ("30M100N28M", 100, 257, "+", 100, 129),
    ("30M100N28M", 100, 257, "-", 100, 129),

    ("3M100N4M",   100, 206, "+", 203, 206),
    ("3M100N4M",   100, 206, "-", 203, 206),

    # clipping should NOT become genomic sequence,
    # but should not prevent correct choice of the longer M-block

    ("2S30M100N28M", 100, 257, "+", 100, 129),
    ("2S30M100N28M", 100, 257, "-", 100, 129),

    ("30M100N28M2S", 100, 257, "+", 100, 129),
    ("30M100N28M2S", 100, 257, "-", 100, 129),

    ("2S3M100N4M", 100, 206, "+", 203, 206),
    ("2S3M100N4M", 100, 206, "-", 203, 206),

    ("3M100N4M2S", 100, 206, "+", 203, 206),
    ("3M100N4M2S", 100, 206, "-", 203, 206),
]

tests += [
    # I не расходует reference и не должен искусственно удлинять genomic span

    # left: 10M + 3I + 10M = 20 aligned bases
    # right: 19M
    # выбираем left
    ("10M3I10M100N19M", 100, 238, "+", 100, 119),
    ("10M3I10M100N19M", 100, 238, "-", 100, 119),

    # D расходует reference, но не увеличивает число aligned query bases
    # left aligned = 20, ref span = 23
    # right aligned = 19
    # выбираем left, coords 100-122
    ("10M3D10M100N19M", 100, 241, "+", 100, 122),
    ("10M3D10M100N19M", 100, 241, "-", 100, 122),

    # = и X считаем как aligned bases и reference-consuming
    # left aligned = 21
    # right aligned = 19
    ("10=1X10=100N19M", 100, 239, "+", 100, 120),
    ("10=1X10=100N19M", 100, 239, "-", 100, 120),

    # справа длиннее
    # left aligned = 19
    # right aligned = 20
    ("19M100N10M2I10M", 100, 238, "+", 219, 238),
    ("19M100N10M2I10M", 100, 238, "-", 219, 238),

    # справа ref span длиннее из-за deletion
    # right aligned = 20, ref span = 23
    ("19M100N10M3D10M", 100, 241, "+", 219, 241),
    ("19M100N10M3D10M", 100, 241, "-", 219, 241),

    # clipping + indel
    ("2S10M3I10M100N19M", 100, 238, "+", 100, 119),
    ("10M3D10M100N19M2S", 100, 241, "-", 100, 122),
]

print(
    f"{'strand':<7} "
    f"{'CIGAR':<18} "
    f"{'coord exp':<14} "
    f"{'coord got':<14} "
    f"RESULT"
)

print("-" * 70)

n_failed = 0

for cigar, start, end, strand, exp_start, exp_end in tests:

    got_start, got_end = coordinate_trimmer_by_cigar(
        EDIT_DIST_TYPE,
        cigar,
        start,
        end,
        strand,
        EXPERIMENT
    )

    expected = (exp_start, exp_end)
    got = (got_start, got_end)

    ok = got == expected

    if not ok:
        n_failed += 1

    print(
        f"{strand:<7} "
        f"{cigar:<18} "
        f"{str(expected):<14} "
        f"{str(got):<14} "
        f"{'OK' if ok else 'FAIL'}"
    )


print()
print(f"Total:  {len(tests)}")
print(f"Passed: {len(tests) - n_failed}")
print(f"Failed: {n_failed}")

if n_failed > 0:
    sys.exit(1)