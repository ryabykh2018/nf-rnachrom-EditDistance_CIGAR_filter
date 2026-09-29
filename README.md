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
1 1 \
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

Passing mates are merged into one DNA interval. RNA output half is filled with `*`

### RNAseq_SE

- only `rna1` is processed
- DNA output half is `*`

### RNAseq_PE

Only `rna1` participates in filtering and final output.

`rna2` is used by the upstream stage of the pipeline for paired-end mapping, proper-pair selection, and unique/multiple mapping determination, but is not used for edit-distance, MAPQ, N filtering, coordinate trimming, or final output.

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

### `filtered_<input_file_name>`

Contains contacts that passed all filters.

Coordinate trimming for a permitted single-`N` CIGAR is applied only to contacts written to this file.

The output is normalized to a fixed RNA–DNA schema:

```text
ATA       -> RNA + DNA
OTA       -> RNA fields = *, DNA populated
RNA-seq   -> RNA populated, DNA fields = *
```

For `RNAseq_PE`, `rna2` is not copied into the DNA half.

### `out_<input_file_name>`

Contains all rejected contacts:

- contacts that passed input validation but failed normal filtering;
- contacts rejected during input validation.

Rejected rows preserve the original upstream coordinates. Coordinate trimming is not applied to rows written to `out_*`.

### `cigar_stat_filtered_<input_file_name>`

Contains the distribution of CIGAR **types** among contacts written to `filtered_*`.

A CIGAR type describes which CIGAR operations are present and how many operation blocks of each type occur. It does **not** store operation lengths.

Examples:

```text
98M                 -> M1
2S96M               -> M1S1
30M100N28M          -> M2N1
10M3D10M100N19M     -> D1M3N1
```

For experiments where both r1 and r2 are evaluated, the two CIGAR types are joined with `-`.

Illustrative ATA statistics file:

```text
CIGAR type (rna-dna)	N	%
M1-M1	                850	85.00
M1S1-M1                100	10.00
M2N1-M1                 50	5.00
```

Here:

- `CIGAR type (...)` is the CIGAR-type class represented by the row;
- `N` is the number of contacts with this CIGAR type;
- `%` is the percentage among all contacts represented in this statistics file.

For downstream modes where only r1 is evaluated (`OTA_SE`, `RNAseq_SE`, `RNAseq_PE`), only the r1-side CIGAR type is reported.

Example:

```text
CIGAR type (rna1)	N	%
M1	                900	90.00
M1S1	                 60	6.00
M2N1	                 40	4.00
```

Accounting invariant:

```text
sum(N in cigar_stat_filtered_*)
=
number of data rows in filtered_*
```

### `cigar_stat_out_<input_file_name>`

Contains the CIGAR-type distribution for contacts that passed input validation but failed **normal filtering**.

Such contacts may fail because of:

- final edit distance above the threshold;
- MAPQ below the threshold;
- an experiment-specific `N` rule;
- the OTA_PE chromosome/distance rule.

This file is a distribution by CIGAR type, not a table of explicit rejection reasons. Therefore, for example, a row with `M1S1` means that contacts with this CIGAR type were rejected, but the CIGAR type itself does not specify which filtering criterion caused rejection.

Illustrative example:

```text
CIGAR type (rna1)	N	%
M1S1	                 70	70.00
M2N1	                 20	20.00
I1M2	                 10	10.00
```

Input-validation failures are intentionally not included in this file.

### `validation_reject_stat_out_<input_file_name>`

Contains counts of contacts rejected because one or more required input fields were invalid or unsupported.

Possible rejection reasons are:

```text
invalid_CIGAR_r1
invalid_CIGAR_r2
invalid_CIGAR_both
invalid_NM
invalid_MAPQ
invalid_coordinates
```

Illustrative example:

```text
validation_reject_reason	N	%
invalid_NM	                12	40.00
invalid_MAPQ	                 8	26.67
invalid_CIGAR_r1	         6	20.00
invalid_coordinates	         4	13.33
```

Here:

- `validation_reject_reason` is the technical validation failure;
- `N` is the number of contacts rejected for that reason;
- `%` is the percentage among all validation rejects.

This separation is intentional:

```text
cigar_stat_out_*
```

describes valid contacts rejected by the normal filtering rules, whereas:

```text
validation_reject_stat_out_*
```

describes malformed or unsupported input.

Rejected-contact accounting invariant:

```text
number of data rows in out_*
=
sum(N in cigar_stat_out_*)
+
sum(N in validation_reject_stat_out_*)
```

### `id_reads_for_ucaRNAs_<input_file_name>`

This file is created only when:

```text
Assembly_of_ucaRNAs = yes
```

and only for ATA experiments:

```text
ATA, not iMARGI
ATA, iMARGI
```

It is not created for OTA or RNA-seq modes.

The file contains two tab-separated columns:

```text
read_id	pairtype
```

Example:

```text
read_id	pairtype
SRR123456.1001	UU
SRR123456.1007	UM
SRR123456.1012	UU
```

The file contains ATA read IDs whose relevant mapped parts passed the downstream filters required for their pairtype.

For `UU` contacts:

```text
UU = uniquely mapped RNA + uniquely mapped DNA
```

both RNA and DNA parts are validated and filtered. The read ID is written to `id_reads_for_ucaRNAs_*` only if both parts pass all applicable checks, including:

```text
valid CIGAR
valid NM
valid MAPQ
valid coordinates
final edit distance <= threshold
MAPQ >= threshold
experiment-specific N rules
```

For ATA, the RNA part may contain at most one `N`, while the DNA part must contain no `N`.

For `UM` contacts:

```text
UM = uniquely mapped RNA + multimapped DNA
```

only the uniquely mapped RNA part is subjected to the downstream CIGAR/NM/MAPQ/coordinate/edit-distance filtering performed by this script. The multimapped DNA part is not treated as a uniquely mapped alignment and is therefore not subjected to the r2 filtering criteria used for `UU`.

A `UM` read ID is written to `id_reads_for_ucaRNAs_*` when the RNA part passes all applicable filters.

Conceptually:

```text
UU -> write read_id if RNA passes AND DNA passes
UM -> write read_id if RNA passes
```

Contacts that fail validation or normal filtering are not included.

The `pairtype` column is retained so downstream ucaRNA assembly can distinguish reads originating from `UU` and `UM` contacts.
## Explorer mode

With:

```text
mode = explorer
```

the script appends six diagnostic CIGAR fields to each output row:

```text
<r1>_cigar_type
<r1>_N_softClipp_bp
<r1>_softClipp_type
<r2>_cigar_type
<r2>_N_softClipp_bp
<r2>_softClipp_type
```

The exact prefixes depend on the experiment:

```text
ATA          -> rna_*  / dna_*
OTA_SE/PE    -> dna1_* / dna2_*
RNAseq_SE/PE -> rna1_* / rna2_*
```

### `*_cigar_type`

A compact description of the operations present in the CIGAR and the number of operation blocks of each type.

Examples:

```text
98M                         -> M1
2S96M                       -> M1S1
10M3I10M                    -> I1M2
30M100N28M                  -> M2N1
1S2M1D2M10N2M1I2M2S        -> D1I1M4N1S2
```

The numbers in `cigar_type` are counts of operation blocks, not numbers of bases. For example:

```text
M2N1
```

means two `M` blocks and one `N` block.

For `ATA, iMARGI`, technical terminal clipping on the original-read 3′ end is removed before the diagnostic CIGAR type and clipping values are calculated.

### `*_N_softClipp_bp`

Contains the total number of terminal clipped bases counted by the script for clipping diagnostics and, when `edit_dist_type = NM + N_softClipp_bp`, for the final edit-distance calculation.

Despite the historical field name, both soft clipping (`S`) and hard clipping (`H`) are counted.

Examples:

```text
98M          -> 0
2S96M        -> 2
5H10S85M     -> 15
2S94M4S      -> 6
```

For double-sided clipping, clipping from the left and right ends is summed in the explorer output.

For `ATA, iMARGI`, clipping at the technical original-read 3′ end is ignored according to strand before this value is calculated.

### `*_softClipp_type`

Describes where the counted terminal clipping is located:

```text
absent  -> no counted terminal S/H
left    -> clipping only at the left end of the CIGAR
right   -> clipping only at the right end of the CIGAR
double  -> clipping at both ends
```

Examples:

```text
98M          -> absent
2S96M        -> left
96M2S        -> right
2S94M4S      -> double
```

### Missing / unused r2 diagnostics

When r2 is not used by the downstream filtering logic, its explorer fields are written as `*`.

This applies to:

- `OTA_SE`;
- `RNAseq_SE`;
- `RNAseq_PE`, where only `rna1` is filtered downstream.

For an ATA `UM` contact, RNA is uniquely mapped and processed normally, while the DNA side is a multimapper. In explorer output the DNA CIGAR type is represented as:

```text
multimapper
```

and the DNA clipping-specific diagnostic fields are `*`.

Explorer mode does not change whether a contact passes or fails filtering. It only appends diagnostic fields to `filtered_*` and `out_*`.
## Assumptions

The script expects the upstream `Bam-to-Contacts` step to provide the required contact-table structure and pairtype classification.

Mapping uniqueness is determined upstream rather than recalculated here. In particular, ATA intentionally accepts both:

```text
UU -> RNA uniquely mapped, DNA uniquely mapped
UM -> RNA uniquely mapped, DNA multimapped
```

For ATA `UU`, both RNA and DNA are subjected to downstream filtering.

For ATA `UM`, only the uniquely mapped RNA part is subjected to downstream filtering; the DNA side remains a multimapper and is not treated as a uniquely mapped r2 alignment.

Where paired-end sequencing/proper-pair selection is applicable, that selection is expected to have been performed upstream. The downstream script does not repeat proper-pair or strand-orientation validation.

The script intentionally does not perform:

```text
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

The final implementation was also regression-checked on representative real datasets for ATA, iMARGI, OTA_SE, OTA_PE, RNAseq_SE and RNAseq_PE.

## Repository checks used before release

```bash
./run_tests.sh
python -m py_compile EditDistance_CIGAR_filter.py
git diff --check
```

A successful full test run ends with all test scripts passing.
