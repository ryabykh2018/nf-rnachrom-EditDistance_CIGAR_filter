# EditDistance_CIGAR_filter

`EditDistance_CIGAR_filter.py` filters raw RNA–chromatin / RNA-seq contact tables after the `Bam-to-Contacts` step and writes a unified RNA–DNA contact table.

This README is the technical reference for running the script: exact input schemas, CLI arguments, output files, validation rules, and tests.

## Requirements

Python dependencies used by the script include:

```text
pandas
matplotlib
```

The input tables are expected to be produced by the corresponding `Bam-to-Contacts` step.

## Supported experiment types

```text
ATA, not iMARGI
ATA, iMARGI
OTA_SE
OTA_PE
RNAseq_SE
RNAseq_PE
```

## Command line

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

The script accepts exactly 12 positional arguments:

| # | Argument | Meaning |
|---:|---|---|
| 1 | `edit_dist_type` | `NM` or `NM + N_softClipp_bp` |
| 2 | `r1_finalThreshold_edit_dist` | maximum final edit distance for r1 |
| 3 | `r2_finalThreshold_edit_dist` | maximum final edit distance for r2 |
| 4 | `r1_mapQ_threshold` | minimum MAPQ for r1 |
| 5 | `r2_mapQ_threshold` | minimum MAPQ for r2 |
| 6 | `OTA_PE_distanceThreshold` | maximum distance between OTA_PE mates |
| 7 | `Assembly_of_ucaRNAs` | `yes` or `no` |
| 8 | `mode` | `not explorer` or `explorer` |
| 9 | `experiment_type` | one of the supported experiment types |
| 10 | `input_file_name` | input table filename |
| 11 | `input_file_path` | input directory |
| 12 | `output_file_path` | output directory |

Numeric thresholds must be non-negative integers.

The script validates CLI values before processing.

## Input requirements

The input is a tab-separated table with a header. A header-only file with zero contacts is supported.

### ATA input header

```text
read_id
ATA_pairtype
rna_chr
rna_start
rna_end
rna_strand
rna_cigar
rna_NM
rna_mapq
dna_chr
dna_start
dna_end
dna_strand
dna_cigar
dna_NM
dna_mapq
rna_secondary_alignments
dna_secondary_alignments
rna_other_tags
dna_other_tags
```

### OTA_SE input header

```text
read_id
OTA_SE_pairtype
dna1_chr
dna1_start
dna1_end
dna1_strand
dna1_cigar
dna1_NM
dna1_mapq
dna2_chr
dna2_start
dna2_end
dna2_strand
dna2_cigar
dna2_NM
dna2_mapq
dna1_secondary_alignments
dna2_secondary_alignments
dna1_other_tags
dna2_other_tags
```

### OTA_PE input header

Same data fields as OTA_SE, with:

```text
OTA_PE_pairtype
```

and both `dna1_*` and `dna2_*` populated.

### RNAseq_SE input header

```text
read_id
RNAseq_SE_pairtype
rna1_chr
rna1_start
rna1_end
rna1_strand
rna1_cigar
rna1_NM
rna1_mapq
rna2_chr
rna2_start
rna2_end
rna2_strand
rna2_cigar
rna2_NM
rna2_mapq
rna1_secondary_alignments
rna2_secondary_alignments
rna1_other_tags
rna2_other_tags
```

### RNAseq_PE input header

Same data fields as RNAseq_SE, with:

```text
RNAseq_PE_pairtype
```

and both RNA mates present in the input.

## Input validation

Before normal filtering, the script validates:

```text
CIGAR
NM
MAPQ
coordinates
```

Supported CIGAR operations:

```text
M I D N S H = X
```

`P` is intentionally unsupported.

Validation-reject categories include:

```text
invalid_CIGAR_r1
invalid_CIGAR_r2
invalid_CIGAR_both
invalid_NM
invalid_MAPQ
invalid_coordinates
```

Validation rejects are written to `out_*` and counted separately from ordinary CIGAR-filter rejects.

## Edit-distance modes

### `NM`

```text
final_edit_distance = NM
```

### `NM + N_softClipp_bp`

```text
final_edit_distance = NM + terminal_clipped_bp
```

Terminal `S` and `H` operations are counted.

For `ATA, iMARGI`, technical clipping at the original-read 3' end is ignored:

```text
strand + -> ignore right terminal clipping
strand - -> ignore left terminal clipping
```

## Experiment-specific `N` rules

| Experiment type | Read | Allowed `N` |
|---|---|---:|
| `ATA, not iMARGI` / `ATA, iMARGI` | RNA | 0 or 1 |
| `ATA, not iMARGI` / `ATA, iMARGI` | DNA | 0 |
| `OTA_SE` | dna1 | 0 |
| `OTA_PE` | dna1 + dna2 | 0 |
| `RNAseq_SE` | rna1 | 0 or 1 |
| `RNAseq_PE` | rna1 only | 0 or 1 |

For exactly one `N`, the selected alignment block is the side with the larger total number of `M`, `=` and `X` bases. Ties select the left block.

Reference span of the selected block includes:

```text
M D = X
```

and excludes:

```text
I S H
```

If there is no `N`, coordinates are not changed.

## Experiment-specific processing

### ATA

- r1 = RNA
- r2 = DNA
- RNA is always processed
- DNA is processed for `UU`
- output contains RNA + DNA

### OTA_SE

- only `dna1` is processed
- any `N` is rejected
- RNA output half is filled with `*`

### OTA_PE

Both mates must pass CIGAR, NM, MAPQ and `N` rules.

Additional conditions:

```text
dna1_chr == dna2_chr
distance <= OTA_PE_distanceThreshold
```

Passing mates are merged into one DNA interval.

### RNAseq_SE

- only `rna1` is processed
- DNA output half is `*`

### RNAseq_PE

Only `rna1` participates in downstream filtering and final output.

`rna2` is upstream pairing / uniqueness evidence and is not used for edit-distance, MAPQ, N filtering, coordinate trimming, or final DNA output.

## Main output schema

Passing contacts are written to a unified RNA–DNA schema:

```text
read_id
pairtype
rna_chr
rna_start
rna_end
rna_strand
rna_cigar
rna_NM
rna_mapq
dna_chr
dna_start
dna_end
dna_strand
dna_cigar
dna_NM
dna_mapq
rna_secondary_alignments
dna_secondary_alignments
rna_other_tags
dna_other_tags
```

Experiment-specific filling:

```text
ATA       -> RNA + DNA
OTA       -> RNA fields = *, DNA populated
RNA-seq   -> RNA populated, DNA fields = *
```

For RNAseq_PE, `rna2` is not copied into the DNA half.

## Output files

```text
filtered_<input_file_name>
out_<input_file_name>
cigar_stat_filtered_<input_file_name>
cigar_stat_out_<input_file_name>
validation_reject_stat_out_<input_file_name>
```

For ATA only, when:

```text
Assembly_of_ucaRNAs = yes
```

the script also writes:

```text
id_reads_for_ucaRNAs_<input_file_name>
```

Diagnostic PNG plots are generated for non-empty NM/clipping/length statistics.

### `filtered_*`

Contains contacts that passed all filters.

Coordinate trimming for single-N contacts is applied only here.

### `out_*`

Contains rejected contacts.

Original input coordinates are preserved.

### Statistics accounting

```text
sum(N in cigar_stat_filtered_*)
=
number of data rows in filtered_*
```

and:

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

six diagnostic fields are appended:

```text
r1_cigar_type
r1_N_softClipp_bp
r1_softClipp_type
r2_cigar_type
r2_N_softClipp_bp
r2_softClipp_type
```

Exact prefixes depend on the experiment header.

Possible clipping types:

```text
absent
left
right
double
```

## Assumptions

The script expects upstream processing to provide the contact-table structure, mapping uniqueness and proper-pair handling where applicable.

The script intentionally does not perform:

```text
CIGAR reference span <-> start/end consistency validation
generic malformed-row column-count validation
OTA_PE strand-orientation validation
```

The output directory is expected to exist.

## Testing

Run the complete test suite:

```bash
./run_tests.sh
```

Current suite:

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

- CIGAR parsing and unsupported/malformed cases;
- soft/hard clipping;
- iMARGI-specific clipping;
- edit-distance modes;
- single-N block selection and coordinates;
- experiment-specific N rules;
- NM/MAPQ/coordinate validation;
- experiment/pairtype matrix;
- RNAseq_PE r1-only behavior;
- OTA_SE / OTA_PE / RNAseq_SE integration;
- ucaRNA output scope;
- statistics accounting;
- header-only input.

The final implementation was also regression-checked on representative real datasets for ATA, iMARGI, OTA_SE, OTA_PE, RNAseq_SE and RNAseq_PE. Observed differences from the old implementation were traced to intended fixes in N handling and clipping-related coordinate logic.

## Repository checks used before release

```bash
./run_tests.sh
python -m py_compile EditDistance_CIGAR_filter.py
git diff --check
```

A successful full test run ends with all test scripts passing.
