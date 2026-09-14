# EditDistance_CIGAR_filter

`EditDistance_CIGAR_filter.py` is a filtering and post-processing step of the `nf-rnachrom` pipeline. It validates mapping-related fields, applies edit-distance, MAPQ and CIGAR-based filters, handles experiment-specific splice-gap rules, adjusts coordinates when a single `N` operation is allowed, and writes a unified RNA–DNA contact table.

## Supported experiment types

- `ATA, not iMARGI`
- `ATA, iMARGI`
- `OTA_SE`
- `OTA_PE`
- `RNAseq_SE`
- `RNAseq_PE`

## Main filtering rules

The script can use either:

```text
NM
```

or:

```text
NM + N_softClipp_bp
```

as the final edit-distance metric.

For `NM + N_softClipp_bp`, terminal `S` and `H` clipping is added to `NM`. For `ATA, iMARGI`, technical clipping at the original-read 3' end is ignored.

Input validation is performed before normal filtering for:

- CIGAR
- NM
- MAPQ
- genomic coordinates

Supported CIGAR operations:

```text
M I D N S H = X
```

`P` is not supported.

### `N` rules

| Experiment type | Read(s) checked | Allowed `N` |
|---|---|---|
| `ATA, not iMARGI` | RNA | 0 or 1 |
| `ATA, not iMARGI` | DNA | 0 |
| `ATA, iMARGI` | RNA | 0 or 1 |
| `ATA, iMARGI` | DNA | 0 |
| `OTA_SE` | dna1 | 0 |
| `OTA_PE` | dna1 + dna2 | 0 |
| `RNAseq_SE` | rna1 | 0 or 1 |
| `RNAseq_PE` | rna1 only | 0 or 1 |

When exactly one `N` is allowed, the script keeps the alignment block with the larger total number of `M`, `=` and `X` bases. Ties are resolved in favor of the left block.

For the selected block, reference span includes:

```text
M D = X
```

and excludes:

```text
I S H N
```

If there is no `N`, genomic coordinates are left unchanged.

## Experiment-specific behavior

### ATA

RNA is always processed. DNA is additionally processed for `UU` contacts.

For `ATA, iMARGI`, technical clipping at the original-read 3' end is excluded from the clipping penalty.

### OTA_SE

Only `dna1` is filtered. The RNA half of the final unified output is filled with `*`.

### OTA_PE

Both DNA mates are filtered independently. Passing contacts must also:

- map to the same chromosome;
- satisfy `distance <= OTA_PE_distanceThreshold`.

The two mates are merged into one DNA interval in the final output.

Proper-pair orientation is expected to be handled upstream.

### RNAseq_SE

Only `rna1` is filtered. The DNA half of the final output is filled with `*`.

### RNAseq_PE

Only `rna1` participates in downstream filtering and output generation.

`rna2` is treated as upstream pairing / uniqueness evidence and is not used for:

- edit-distance filtering;
- MAPQ filtering;
- `N` filtering;
- coordinate trimming;
- final RNA–DNA output.

## Input assumptions

The input is a tab-separated contacts table with a header.

A header-only file with zero contacts is supported.

The script assumes that upstream pipeline stages have already handled basic table construction and, where applicable, proper-pair and uniqueness logic.

The script intentionally does not validate:

- CIGAR reference span against `start/end`;
- total column count of malformed rows;
- OTA_PE strand orientation.

## Usage

```bash
python EditDistance_CIGAR_filter.py \
"NM + N_softClipp_bp" \
2 2 \
0 0 \
200 \
"yes" \
"not explorer" \
"ATA, not iMARGI" \
"input_contacts.tab.rc" \
"/path/to/input/" \
"/path/to/output/"
```

The script accepts 12 positional arguments:

| # | Argument | Description |
|---:|---|---|
| 1 | `edit_dist_type` | `NM` or `NM + N_softClipp_bp` |
| 2 | `r1_finalThreshold_edit_dist` | Maximum final edit distance for r1 |
| 3 | `r2_finalThreshold_edit_dist` | Maximum final edit distance for r2 |
| 4 | `r1_mapQ_threshold` | Minimum MAPQ for r1 |
| 5 | `r2_mapQ_threshold` | Minimum MAPQ for r2 |
| 6 | `OTA_PE_distanceThreshold` | Maximum allowed distance between OTA_PE mates |
| 7 | `Assembly_of_ucaRNAs` | `yes` or `no` |
| 8 | `mode` | `not explorer` or `explorer` |
| 9 | `experiment_type` | One of the supported experiment types |
| 10 | `input_file_name` | Input contacts-table filename |
| 11 | `input_file_path` | Input directory |
| 12 | `output_file_path` | Output directory |

All numeric thresholds must be non-negative integers.

Input and output directory arguments are expected to end with `/`.

## Output files

The main outputs are:

```text
filtered_<input_file_name>
out_<input_file_name>
cigar_stat_filtered_<input_file_name>
cigar_stat_out_<input_file_name>
validation_reject_stat_out_<input_file_name>
```

For ATA experiments, when:

```text
Assembly_of_ucaRNAs = yes
```

the script also creates:

```text
id_reads_for_ucaRNAs_<input_file_name>
```

Diagnostic PNG plots are generated from edit-distance, clipping and length statistics when the corresponding statistics are non-empty.

### `filtered_*`

Contains contacts that passed all filters.

The output is normalized to a fixed RNA–DNA schema:

- ATA: RNA + DNA populated;
- OTA: RNA half is `*`, DNA populated;
- RNA-seq: RNA populated, DNA half is `*`.

### `out_*`

Contains rejected contacts.

Rejected rows preserve their original input coordinates. Coordinate trimming is applied only to contacts written to `filtered_*`.

### `cigar_stat_filtered_*`

CIGAR-type distribution for contacts that passed filtering.

Accounting invariant:

```text
sum(N in cigar_stat_filtered_*) == number of data rows in filtered_*
```

### `cigar_stat_out_*`

CIGAR-type distribution for contacts that passed input validation but failed normal filtering.

Validation failures are not written here.

### `validation_reject_stat_out_*`

Counts input-validation failures such as:

```text
invalid_CIGAR_r1
invalid_CIGAR_r2
invalid_CIGAR_both
invalid_NM
invalid_MAPQ
invalid_coordinates
```

Accounting invariant:

```text
number of data rows in out_*
=
sum(N in cigar_stat_out_*)
+
sum(N in validation_reject_stat_out_*)
```

## Explorer mode

With:

```text
mode = explorer
```

additional CIGAR diagnostics are appended to output rows, including:

- CIGAR type;
- number of clipped bases;
- clipping type.

Possible clipping types are:

```text
absent
left
right
double
```

## Tests

Run the full test suite with:

```bash
./run_tests.sh
```

Current tests:

```text
test_cigar_imargi.py
test_cigar_single_N.py
test_cigar_N_validation.py
test_N_integration.py
test_edit_distance_clipping.py
test_cigar_type.py
test_edit_distance_integration.py
test_single_N_coordinates_integration.py
test_imargi_single_N_integration.py
test_cigar_parser_edges.py
test_invalid_cigar_integration.py
test_invalid_NM_integration.py
test_invalid_MAPQ_integration.py
test_invalid_coordinates_integration.py
test_experiment_pairtype_matrix.py
test_rnaseq_pe_r1_only.py
test_ota_pe_filtering.py
test_ota_se_filtering.py
test_rnaseq_se_filtering.py
test_ucaRNA_output_scope.py
test_statistics_accounting.py
test_header_only_input.py
```

The suite covers:

- CIGAR parsing and edge cases;
- soft/hard clipping;
- iMARGI clipping rules;
- edit-distance calculation;
- single-`N` coordinate trimming;
- experiment-specific `N` rules;
- malformed CIGAR / NM / MAPQ / coordinates;
- ATA / OTA / RNA-seq experiment logic;
- RNAseq_PE r1-only behavior;
- ucaRNA output scope;
- statistics accounting;
- header-only input.

## Documentation

A detailed Russian-language description of the algorithm, experiment-specific rules, output files and examples is provided separately in:

```text
EditDistance_CIGAR_filter_documentation_RU.docx
```
