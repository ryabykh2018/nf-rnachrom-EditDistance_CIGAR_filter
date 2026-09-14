import pandas as pd
from matplotlib import pyplot as plt
from collections import Counter
import sys
import re


def headerParser(experiment_type):
    SRR_ID = 'read_id'
    if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
        pairtype = 'ATA_pairtype'
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'rna_chr', 'rna_start', 'rna_end', 'rna_strand', 'rna_cigar', 'rna_secondary_alignments', 'rna_other_tags', 'rna_NM', 'rna_mapq'
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'dna_chr', 'dna_start', 'dna_end', 'dna_strand', 'dna_cigar', 'dna_secondary_alignments', 'dna_other_tags', 'dna_NM', 'dna_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "rna_cigar_type", "rna_N_softClipp_bp", "rna_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "dna_cigar_type", "dna_N_softClipp_bp", "dna_softClipp_type"

    elif experiment_type in ['OTA_SE', 'OTA_PE']:
        pairtype = experiment_type + '_pairtype'
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'dna1_chr', 'dna1_start', 'dna1_end', 'dna1_strand', 'dna1_cigar', 'dna1_secondary_alignments', 'dna1_other_tags', 'dna1_NM', 'dna1_mapq'
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'dna2_chr', 'dna2_start', 'dna2_end', 'dna2_strand', 'dna2_cigar', 'dna2_secondary_alignments', 'dna2_other_tags', 'dna2_NM', 'dna2_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "dna1_cigar_type", "dna1_N_softClipp_bp", "dna1_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "dna2_cigar_type", "dna2_N_softClipp_bp", "dna2_softClipp_type"

    elif experiment_type in ['RNAseq_SE', 'RNAseq_PE']:
        pairtype = experiment_type + '_pairtype'
        r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r1_mapQ = 'rna1_chr', 'rna1_start', 'rna1_end', 'rna1_strand', 'rna1_cigar', 'rna1_secondary_alignments', 'rna1_other_tags', 'rna1_NM', 'rna1_mapq'
        r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r2_mapQ = 'rna2_chr', 'rna2_start', 'rna2_end', 'rna2_strand', 'rna2_cigar', 'rna2_secondary_alignments', 'rna2_other_tags', 'rna2_NM', 'rna2_mapq'

        r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = "rna1_cigar_type", "rna1_N_softClipp_bp", "rna1_softClipp_type"
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = "rna2_cigar_type", "rna2_N_softClipp_bp", "rna2_softClipp_type"

    else:
        raise ValueError(f"Unsupported experiment type: {experiment_type}")

    return SRR_ID, pairtype, r1_chr, r1_start, r1_end, r1_strand, r1_cigar, r1_secondary_alignments, r1_other_tags, r1_NM, r2_chr, r2_start, r2_end, r2_strand, r2_cigar, r2_secondary_alignments, r2_other_tags, r2_NM, r1_mapQ, r2_mapQ, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type #r1_final_edit_dist, , r2_final_edit_dist




def CIGAR_field_classifier_parent(pairtype, edit_dist_type, r1_cigar, r1_strand, r1_NM, r2_cigar, r2_strand, r2_NM, experiment_type):
    r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type = CIGAR_field_classifier(r1_cigar, experiment_type, r1_strand)
    r1_final_edit_dist = calculate_edit_distance(r1_NM, r1_N_softClipp_bp, edit_dist_type)
    # r2 is filtered downstream only for ATA UU contacts and OTA_PE.
    # For RNAseq_PE, r2 is upstream pairing/uniqueness evidence only; final filtering and coordinates are based on r1.
    if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']):
        r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type = CIGAR_field_classifier(r2_cigar, experiment_type, r2_strand)
        r2_final_edit_dist = calculate_edit_distance(r2_NM, r2_N_softClipp_bp, edit_dist_type)
    else:
        r2_cigar_type = "multimapper" if pairtype == 'UM' else "*"
        r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = "*", "*", "*"
    return r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist



def CIGAR_field_classifier(CIGAR_field, experiment_type, strand):
    CLIPPING_SYMBOLS = {'S', 'H'}

    # Parse CIGAR into [(length, operation), ...]
    tokens = parse_cigar_tokens(CIGAR_field)

    if tokens is None:
        raise ValueError(f"Invalid or unsupported CIGAR: {CIGAR_field}")

    # For ATA iMARGI, ignore technical clipping at the original read 3' end:
    # strand "+" -> original-read 3' end is the right CIGAR end
    # strand "-" -> original-read 3' end is the left CIGAR end
    if experiment_type == "ATA, iMARGI":
        if strand == "+":
            while tokens and tokens[-1][1] in CLIPPING_SYMBOLS:
                tokens.pop()

        elif strand == "-":
            while tokens and tokens[0][1] in CLIPPING_SYMBOLS:
                tokens.pop(0)

    # Count clipping remaining after the iMARGI-specific removal
    left_clip_bp = 0
    right_clip_bp = 0

    i = 0
    while i < len(tokens) and tokens[i][1] in CLIPPING_SYMBOLS:
        left_clip_bp += tokens[i][0]
        i += 1

    i = len(tokens) - 1
    while i >= 0 and tokens[i][1] in CLIPPING_SYMBOLS:
        right_clip_bp += tokens[i][0]
        i -= 1

    if left_clip_bp and right_clip_bp:
        softClipp_type = "double"
        N_softClipp_bp = [left_clip_bp, right_clip_bp]

    elif left_clip_bp:
        softClipp_type = "left"
        N_softClipp_bp = [left_clip_bp]

    elif right_clip_bp:
        softClipp_type = "right"
        N_softClipp_bp = [right_clip_bp]

    else:
        softClipp_type = "absent"
        N_softClipp_bp = [0]

    # Preserve previous cigar_type idea:
    # e.g. 1S2M1D2M10N2M1I2M2S -> D1I1M4N1S2
    cigar_ops = [op for _, op in tokens]

    counts = Counter(cigar_ops)
    cigar_type = ''.join(
        op + str(counts[op])
        for op in sorted(counts)
    )

    return cigar_type, N_softClipp_bp, softClipp_type



def check_distance_and_chromosomes(r1_chr, r1_start, r1_end, r2_chr, r2_start, r2_end, OTA_PE_distanceThreshold):
    distance = max(int(r1_start) - int(r2_end) - 1, int(r2_start) - int(r1_end) - 1)
    return (distance <= OTA_PE_distanceThreshold) and (r1_chr == r2_chr)




def parse_cigar_tokens(cigar):
    if not cigar or cigar == "*":
        return None

    tokens = re.findall(r'(\d+)([MIDNSH=X])', cigar)

    # Make sure the parsed tokens cover the entire CIGAR string
    reconstructed = ''.join(length + op for length, op in tokens)

    if reconstructed != cigar:
        return None

    # Zero-length CIGAR operations are invalid
    if any(int(length) <= 0 for length, op in tokens):
        return None

    return [(int(length), op) for length, op in tokens]



def update_counter(counter, key):
    """Increment a dictionary counter by key."""
    counter[key] = counter.get(key, 0) + 1


def cigar_explorer_fields(
    r1_cigar_type,
    r1_N_softClipp_bp,
    r1_softClipp_type,
    r2_cigar_type="*",
    r2_N_softClipp_bp="*",
    r2_softClipp_type="*"
):
    """Build the six CIGAR diagnostic fields used in explorer mode."""

    r1_clipped_bp = str(sum(r1_N_softClipp_bp))

    if r2_N_softClipp_bp == "*":
        r2_clipped_bp = "*"
    else:
        r2_clipped_bp = str(sum(r2_N_softClipp_bp))

    return [
        r1_cigar_type,
        r1_clipped_bp,
        r1_softClipp_type,
        r2_cigar_type,
        r2_clipped_bp,
        r2_softClipp_type
    ]


def write_validation_reject(
    output_file,
    line,
    mode,
    explorer_fields,
    validation_statistics_dict,
    reason
):
    """
    Write a contact rejected during input validation
    and update validation-rejection statistics.
    """

    if mode == "explorer":
        output = line + "\t" + "\t".join(map(str, explorer_fields))
    else:
        output = line

    output_file.write(output + "\n")
    update_counter(validation_statistics_dict, reason)



def validate_NM(NM):
    try:
        value = float(NM)
    except (TypeError, ValueError):
        return False

    # NM must be a non-negative integer value
    if not value.is_integer():
        return False

    if value < 0:
        return False

    return True



def validate_MAPQ(mapq):
    try:
        value = int(mapq)
    except (TypeError, ValueError):
        return False

    if value < 0:
        return False

    return True



def validate_coordinates(start, end):
    try:
        start_value = int(start)
        end_value = int(end)
    except (TypeError, ValueError):
        return False

    if start_value <= 0:
        return False

    if end_value <= 0:
        return False

    if start_value > end_value:
        return False

    return True




def calculate_edit_distance(NM, N_softClipp_bp, edit_dist_type):
    if edit_dist_type == "NM + N_softClipp_bp":
        return int(float(NM)) + sum(N_softClipp_bp)
    elif edit_dist_type == "NM":
        return int(float(NM))
    else:
        raise ValueError(f"Unsupported edit distance type: {edit_dist_type}")



def validate_cigar_N(cigar, experiment_type, read_role):
    n_count = cigar.count("N")

    if experiment_type in ["ATA, not iMARGI", "ATA, iMARGI"]:
        if read_role == "rna":
            return n_count <= 1
        if read_role == "dna":
            return n_count == 0
        return False

    if experiment_type in ["OTA_SE", "OTA_PE"]:
        return n_count == 0

    if experiment_type in ["RNAseq_SE", "RNAseq_PE"]:
        return n_count <= 1

    return False



def coordinate_trimmer_by_cigar(edit_dist_type, cigar, start, end,  strand, experiment_type):

    # No splice gap -> pysam coordinates already describe the actually aligned genomic interval.
    if "N" not in cigar:
        return start, end

    tokens = parse_cigar_tokens(cigar)

    if tokens is None:
        raise ValueError(f"Invalid or unsupported CIGAR: {cigar}")

    # For now we support exactly one N
    n_positions = [i for i, (_, op) in enumerate(tokens) if op == "N"]

    if len(n_positions) != 1:
        # multi-N should be filtered separately from now on
        return start, end

    n_idx = n_positions[0]

    left_tokens = tokens[:n_idx]
    right_tokens = tokens[n_idx + 1:]

    # How many actual aligned query bases are in each block.
    # "I" is not counted as matched/aligned-to-reference bases,
    # "S/H" is also not counted.
    aligned_ops = {"M", "=", "X"}

    left_aligned_len = sum(
        length for length, op in left_tokens if op in aligned_ops
    )

    right_aligned_len = sum(
        length for length, op in right_tokens if op in aligned_ops
    )

    # Reference span of each block.
    # D occupies reference coordinates.
    # I/S/H do not occupy reference coordinates.
    reference_ops = {"M", "D", "=", "X"}

    left_ref_len = sum(
        length for length, op in left_tokens if op in reference_ops
    )

    right_ref_len = sum(
        length for length, op in right_tokens if op in reference_ops
    )

    # pysam start/end include:
    # left block + N + right block
    #
    # Select the block with the largest number of aligned bases.
    # If there are ties, keep the old logic: select the left one.
    if left_aligned_len >= right_aligned_len:
        start_new = start
        end_new = start + left_ref_len - 1

    else:
        start_new = end - right_ref_len + 1
        end_new = end

    return start_new, end_new





def output_string_modifier(output_string, mode, header_dict, edit_dist_type, experiment_type, softClipp_NM_statistics_dict, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type):

    contact = output_string.split('\t')

    SRR_ID                  = contact[header_dict[SRR_ID_header]]#
    pairtype                = contact[header_dict[pairtype_header]]#
    r1_chr                  = contact[header_dict[r1_chr_header]]
    r1_start                = contact[header_dict[r1_start_header]]#
    r1_end                  = contact[header_dict[r1_end_header]]#
    r1_strand               = contact[header_dict[r1_strand_header]]#
    r1_cigar                = contact[header_dict[r1_cigar_header]]#
    r1_NM                   = contact[header_dict[r1_NM_header]]
    r1_mapQ                 = contact[header_dict[r1_mapQ_header]]
    r2_chr                  = contact[header_dict[r2_chr_header]]
    r2_start                = contact[header_dict[r2_start_header]]#
    r2_end                  = contact[header_dict[r2_end_header]]#
    r2_strand               = contact[header_dict[r2_strand_header]]#
    r2_cigar                = contact[header_dict[r2_cigar_header]]#
    r2_NM                   = contact[header_dict[r2_NM_header]]
    r2_mapQ                 = contact[header_dict[r2_mapQ_header]]
    r1_secondary_alignments = contact[header_dict[r1_secondary_alignments_header]]
    r2_secondary_alignments = contact[header_dict[r2_secondary_alignments_header]]
    r1_other_tags           = contact[header_dict[r1_other_tags_header]]
    r2_other_tags           = contact[header_dict[r2_other_tags_header]]

    r1_start_new, r1_end_new = coordinate_trimmer_by_cigar(edit_dist_type, r1_cigar, int(r1_start), int(r1_end), r1_strand, experiment_type)
    if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']):
        r2_start_new, r2_end_new = coordinate_trimmer_by_cigar(edit_dist_type, r2_cigar, int(r2_start), int(r2_end), r2_strand, experiment_type)
    else:
        r2_start_new, r2_end_new = r2_start, r2_end

    r1, r2 = r1_chr_header.split("_")[0], r2_chr_header.split("_")[0]
    if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
        output = "\t".join([
            SRR_ID, pairtype,
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            r2_chr, str(r2_start_new), str(r2_end_new), r2_strand, r2_cigar, r2_NM, r2_mapQ,
            r1_secondary_alignments, r2_secondary_alignments, r1_other_tags, r2_other_tags
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)
        if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']):
            length_statistics(softClipp_NM_statistics_dict, r2_start_new, r2_end_new, r2)

    elif experiment_type in ['RNAseq_SE', 'RNAseq_PE']:
        output = "\t".join([
            SRR_ID, pairtype,
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            "*", "*", "*", "*", "*", "*", "*",
            r1_secondary_alignments, "*", r1_other_tags, "*"
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)

    elif experiment_type == 'OTA_SE':
        output = "\t".join([
            SRR_ID, pairtype,
            "*", "*", "*", "*", "*", "*", "*",
            r1_chr, str(r1_start_new), str(r1_end_new), r1_strand, r1_cigar, r1_NM, r1_mapQ,
            "*", r1_secondary_alignments, "*", r1_other_tags
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, r1)

    else: #experiment_type == 'OTA_PE'
        r1_start_new = min(r1_start_new, r2_start_new)
        r1_end_new = max(r1_end_new, r2_end_new)
        output = "\t".join([
            SRR_ID, pairtype,
            "*", "*", "*", "*", "*", "*", "*",
            r1_chr, str(r1_start_new), str(r1_end_new), ";".join(["dna1:" + r1_strand, "dna2:" + r2_strand]), ";".join(["dna1:" + r1_cigar, "dna2:" + r2_cigar]), ";".join(["dna1:" + r1_NM, "dna2:" + r2_NM]), ";".join(["dna1:" + r1_mapQ,"dna2:" + r2_mapQ]),
            "*", "*", ";".join(["dna1:" + r1_secondary_alignments, "dna2:" + r2_secondary_alignments]),  ";".join(["dna1:" + r1_other_tags, "dna2:" + r2_other_tags])
        ])
        length_statistics(softClipp_NM_statistics_dict, r1_start_new, r1_end_new, 'dna')

    if mode == "explorer":
        if r2_N_softClipp_bp != "*":
            r2_N_softClipp_bp = str(sum(r2_N_softClipp_bp))
        output = output + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type])
    return output


def length_statistics(softClipp_NM_statistics_dict, start, end, readType):
    length = end - start + 1
    if length in softClipp_NM_statistics_dict["{} length".format(readType)].keys():
        softClipp_NM_statistics_dict["{} length".format(readType)][length] += 1
    else:
        softClipp_NM_statistics_dict["{} length".format(readType)][length] = 1




def softClipp_NM_statistics(softClipp_NM_statistics_dict, r1_N_softClipp_bp, r2_N_softClipp_bp, r1_NM, r2_NM, r1_chr_header, r2_chr_header):
    r1 = r1_chr_header.split("_")[0]
    r2 = r2_chr_header.split("_")[0]
    if r2_N_softClipp_bp != "*":
        N_softClipp_type = str(sum(r1_N_softClipp_bp)) + "-" + str(sum(r2_N_softClipp_bp))
        if N_softClipp_type in softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)].keys():
            softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)][N_softClipp_type] += 1
        else:
            softClipp_NM_statistics_dict["N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2)][N_softClipp_type] = 1

        NM_NM_type = str(r1_NM) + "-" + str(r2_NM)
        if NM_NM_type in softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)].keys():
            softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)][NM_NM_type] += 1
        else:
            softClipp_NM_statistics_dict["NM({r1}-{r2})".format(r1=r1, r2=r2)][NM_NM_type] = 1
    else:
        NM_softClipp_type = str(r1_NM) + "-" + str(sum(r1_N_softClipp_bp))
        if NM_softClipp_type in softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)].keys():
            softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)][NM_softClipp_type] += 1
        else:
            softClipp_NM_statistics_dict["NM-N_softClipp_bp({r1})".format(r1=r1)][NM_softClipp_type] = 1


def CIGAR_statistics(CIGAR_statistics_dict, experiment_type, r1_cigar_type, r2_cigar_type):
    if experiment_type in ['OTA_SE', 'RNAseq_SE', 'RNAseq_PE']:
        cigar_type = r1_cigar_type
    else:
        cigar_type = r1_cigar_type + "-" + r2_cigar_type
    if cigar_type in CIGAR_statistics_dict.keys():
        CIGAR_statistics_dict[cigar_type] += 1
    else:
        CIGAR_statistics_dict[cigar_type] = 1


def save_cigar_statistics(CIGAR_statistics_dict, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, filtered_or_filtered_out):
    if experiment_type in ['OTA_SE', 'RNAseq_SE', 'RNAseq_PE']:
        CIGAR_statistics_df = pd.DataFrame.from_dict({'CIGAR type ({r1})'.format(r1=r1_chr_header.split('_')[0]): list(CIGAR_statistics_dict.keys()), 'N': list(CIGAR_statistics_dict.values())})
    else:
        CIGAR_statistics_df = pd.DataFrame.from_dict({'CIGAR type ({r1}-{r2})'.format(r1=r1_chr_header.split('_')[0], r2=r2_chr_header.split('_')[0]): list(CIGAR_statistics_dict.keys()), 'N': list(CIGAR_statistics_dict.values())})
    CIGAR_statistics_df['%'] = (100 * CIGAR_statistics_df['N'] / CIGAR_statistics_df['N'].sum()).round(2)
    CIGAR_statistics_df.sort_values(by='N', ascending=False).to_csv(output_file_path + 'cigar_stat_{}_'.format(filtered_or_filtered_out) + input_file_name, index=False, sep='\t')


def save_validation_reject_statistics(validation_statistics_dict, output_file_path, input_file_name):
    """Save counts of contacts rejected during input validation."""
    validation_statistics_df = pd.DataFrame.from_dict({"validation_reject_reason": list(validation_statistics_dict.keys()), "N": list(validation_statistics_dict.values())})
    if not validation_statistics_df.empty:
        total = validation_statistics_df["N"].sum()
        validation_statistics_df["%"] = (100 * validation_statistics_df["N"] / total).round(2)
        validation_statistics_df = validation_statistics_df.sort_values(by="N",ascending=False)
    else:
        validation_statistics_df["%"] = pd.Series(dtype=float)
    validation_statistics_df.to_csv( output_file_path + "validation_reject_stat_out_" + input_file_name, index=False, sep="\t" )




def plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_dict, experiment_type, r1, r2, dataType, path):
    n = 0
    for i in ['N_softClipp_bp({r1}-{r2})'.format(r1=r1, r2=r2), 'NM({r1}-{r2})'.format(r1=r1, r2=r2), 'NM-N_softClipp_bp({})'.format(r1)]:
        subdict = N_softClipp_bp_or_NM_dict[i]
        if len(subdict.keys()) != 0:
            n += 1
            s = {i: [], 'N': []}
            for key in subdict:
                s[i].append(key)
                s['N'].append(subdict[key])
            df = pd.DataFrame.from_dict(s)
            df['{}'.format(r1)] = df[i].apply(lambda x: int(x.split("-")[0]))
            df['{}'.format(r2)] = df[i].apply(lambda x: int(x.split("-")[1]))

            r1_N_softClipp_or_NM = df[['{}'.format(r1), 'N']].groupby(['{}'.format(r1)]).sum().reset_index()
            r2_N_softClipp_or_NM = df[['{}'.format(r2), 'N']].groupby(['{}'.format(r2)]).sum().reset_index()

            min_r1_N_softClipp_or_NM = r1_N_softClipp_or_NM['{}'.format(r1)].min()
            max_r1_N_softClipp_or_NM = r1_N_softClipp_or_NM['{}'.format(r1)].max()
            min_r2_N_softClipp_or_NM = r2_N_softClipp_or_NM['{}'.format(r2)].min()
            max_r2_N_softClipp_or_NM = r2_N_softClipp_or_NM['{}'.format(r2)].max()

            fig = plt.figure()
            fig.set_figheight(8)
            fig.set_figwidth(8)
            ax1 = plt.subplot2grid(shape=(3, 3), loc=(0, 0), colspan=2)
            ax2 = plt.subplot2grid(shape=(3, 3), loc=(1, 0), colspan=2, rowspan=2)
            ax3 = plt.subplot2grid(shape=(3, 3), loc=(1, 2), rowspan=2)


            ax1.bar(r1_N_softClipp_or_NM['{}'.format(r1)].values, r1_N_softClipp_or_NM['N'].values)
            ax1.set_xlim([min_r1_N_softClipp_or_NM - 1, max_r1_N_softClipp_or_NM + 1])
            if r1_N_softClipp_or_NM['N'].max() > 100:
                ax1.set_yscale('log')
            ax1.set_ylabel('N')
            ax1.tick_params(axis='x', which='both', labelcolor="w")


            ax2.scatter(df['{}'.format(r1)].values, df['{}'.format(r2)].values, s=50)
            ax2.set_xlim([min_r1_N_softClipp_or_NM - 1, max_r1_N_softClipp_or_NM + 1])
            ax2.set_ylim([min_r2_N_softClipp_or_NM - 1, max_r2_N_softClipp_or_NM + 1])
            if i == 'NM-N_softClipp_bp({})'.format(r1):
                pairtype = " (pairtype = UM)" if experiment_type in ['ATA, iMARGI', 'ATA, not iMARGI'] else ""
                ax2.set_xlabel('{}_NM'.format(r1))
                ax2.set_ylabel('{}_N_softClipp_bp'.format(r1))
                fig.suptitle("Edit distance and soft-clipped bases in {r1}{pairtype} ({dataType})".format(r1=r1, pairtype=pairtype, dataType=dataType))
            else:
                pairtype = " (pairtype = UU)" if experiment_type in ['ATA, iMARGI', 'ATA, not iMARGI'] else ""
                if i.split("(")[0] == 'N_softClipp_bp':
                    j = "Soft-clipped bases in {r1} and {r2}{pairtype} ({dataType})".format(r1=r1, r2=r2, pairtype=pairtype, dataType=dataType)
                else:
                    j = "Edit distance in {r1} and {r2}{pairtype} ({dataType})".format(r1=r1, r2=r2, pairtype=pairtype, dataType=dataType)
                ax2.set_xlabel('{r1}_{i}'.format(r1=r1, i=i.split("(")[0]))
                ax2.set_ylabel('{r2}_{i}'.format(r2=r2, i=i.split("(")[0]))
                fig.suptitle(j)


            ax3.barh(r2_N_softClipp_or_NM['{}'.format(r2)].values, r2_N_softClipp_or_NM['N'].values)
            ax3.set_ylim([min_r2_N_softClipp_or_NM - 1, max_r2_N_softClipp_or_NM + 1])
            if r2_N_softClipp_or_NM['N'].max() > 100:
                ax3.set_xscale('log')
            ax3.set_xlabel('N')
            ax3.tick_params(axis='y', which='both', labelcolor="w")

            fig.savefig('{path}{dataType}_{n}.png'.format(path=path, dataType=dataType, n=n), dpi=360)
            plt.close(fig)
    if dataType == "filtered_data":
        r1 = r1 if experiment_type != "OTA_PE" else "dna"
        for i in ['{} length'.format(r1), '{} length'.format(r2)]:
            if len(N_softClipp_bp_or_NM_dict[i].keys()) != 0:
                g = plt.figure()
                plt.hist(N_softClipp_bp_or_NM_dict[i].keys(), bins=50, weights=N_softClipp_bp_or_NM_dict[i].values())
                plt.yscale('log')
                g.suptitle('Distribution of {} part lengths ({dataType})'.format(i.split(" ")[0], dataType=dataType))
                plt.xlabel(i)
                plt.ylabel('N')
                g.savefig('{path}{dataType}_{i}_length.png'.format(path=path, dataType=dataType, i=i.split(" ")[0]), dpi=360)
                plt.close(g)






def editDistance_and_CIGAR_filter(edit_dist_type, r1_finalThreshold_edit_dist, r2_finalThreshold_edit_dist, r1_mapQ_threshold, r2_mapQ_threshold, OTA_PE_distanceThreshold, Assembly_of_ucaRNAs, mode, experiment_type, input_file_name, input_file_path, output_file_path):

    ucaRNA_enabled = (Assembly_of_ucaRNAs == "yes" and experiment_type in ["ATA, not iMARGI", "ATA, iMARGI"])

    raw_contacts = open(input_file_path + input_file_name, 'r')
    CIGAR_filtered_contacts = open(output_file_path + 'filtered_' + input_file_name, 'w')
    CIGAR_filtered_out_contacts = open(output_file_path + 'out_' + input_file_name, 'w')
    if ucaRNA_enabled:
        id_reads_for_ucaRNAs = open(output_file_path + 'id_reads_for_ucaRNAs_' + input_file_name, 'w')

    CIGAR_statistics_for_filtered_data = {}
    CIGAR_statistics_for_filtered_out_data = {}
    validation_reject_statistics = {}

    header_dict = {}
    count = 0
    for line in raw_contacts:
        count += 1
        line = line.strip()
        if count == 1:
            SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header = headerParser(experiment_type) #, mode
            if mode != "explorer":
                header_output_filtered_out = line
                header_output_filtered = line.replace("dna1","rna").replace("dna2","dna").replace("rna1","rna").replace("rna2","dna").replace("ATA_pairtype","pairtype").replace("OTA_SE_pairtype","pairtype").replace("OTA_PE_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype")
            else:
                header_output_filtered_out = line + "\t" + "\t".join([r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header])
                header_output_filtered = line.replace("dna1","rna").replace("dna2","dna").replace("rna1","rna").replace("rna2","dna").replace("ATA_pairtype","pairtype").replace("OTA_SE_pairtype","pairtype").replace("OTA_PE_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype") + "\t" + "\t".join([r1_cigar_type_header, r1_N_softClipp_bp_header, r1_softClipp_type_header, r2_cigar_type_header, r2_N_softClipp_bp_header, r2_softClipp_type_header])
            r1, r2 = r1_chr_header.split("_")[0], r2_chr_header.split("_")[0]
            r1_final = "dna" if experiment_type == 'OTA_PE' else r1
            N_softClipp_bp_or_NM_for_filtered_data = {"N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM-N_softClipp_bp({})".format(r1): {},
                                                      "{} length".format(r1_final): {},
                                                      "{} length".format(r2): {}}

            N_softClipp_bp_or_NM_for_filtered_out_data = {"N_softClipp_bp({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM({r1}-{r2})".format(r1=r1, r2=r2): {}, "NM-N_softClipp_bp({})".format(r1): {}}
            CIGAR_filtered_contacts.write(header_output_filtered + "\n")
            CIGAR_filtered_out_contacts.write(header_output_filtered_out + "\n")
            header = line.split('\t')
            for i in range(len(header)):
                header_dict[header[i]] = i
            if ucaRNA_enabled:
                id_reads_for_ucaRNAs.write(SRR_ID_header + "\t" + pairtype_header.replace("ATA_pairtype","pairtype").replace("RNAseq_SE_pairtype","pairtype").replace("RNAseq_PE_pairtype","pairtype") + "\n")
        else:
            contact = line.split('\t')
            pairtype = contact[header_dict[pairtype_header]]

            if (pairtype == 'UU') and (experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI', 'OTA_PE']):
                # Validate both CIGAR strings before classification
                r1_cigar_raw = contact[header_dict[r1_cigar_header]]
                r2_cigar_raw = contact[header_dict[r2_cigar_header]]

                r1_cigar_valid = parse_cigar_tokens(r1_cigar_raw) is not None
                r2_cigar_valid = parse_cigar_tokens(r2_cigar_raw) is not None

                if not (r1_cigar_valid and r2_cigar_valid):
                    if r1_cigar_valid:
                        r1_type, r1_clip, r1_clip_type = CIGAR_field_classifier( r1_cigar_raw, experiment_type, contact[header_dict[r1_strand_header]] )
                        r1_extra = [ r1_type, str(sum(r1_clip)), r1_clip_type ]
                    else:
                        r1_extra = [ "invalid", "*", "*" ]

                    if r2_cigar_valid:
                        r2_type, r2_clip, r2_clip_type = CIGAR_field_classifier( r2_cigar_raw, experiment_type, contact[header_dict[r2_strand_header]] )
                        r2_extra = [ r2_type, str(sum(r2_clip)), r2_clip_type ]
                    else:
                        r2_extra = ["invalid","*","*"]

                    if not r1_cigar_valid and not r2_cigar_valid:
                        invalid_cigar_type = "invalid_CIGAR_both"
                    elif not r1_cigar_valid:
                        invalid_cigar_type = "invalid_CIGAR_r1"
                    else:
                        invalid_cigar_type = "invalid_CIGAR_r2"

                    write_validation_reject( CIGAR_filtered_out_contacts, line, mode, r1_extra + r2_extra, validation_reject_statistics, invalid_cigar_type )
                    continue

                # Validate NM values before edit-distance calculation
                r1_NM_raw = contact[header_dict[r1_NM_header]]
                r2_NM_raw = contact[header_dict[r2_NM_header]]
                r1_NM_valid = validate_NM(r1_NM_raw)
                r2_NM_valid = validate_NM(r2_NM_raw)
                if not (r1_NM_valid and r2_NM_valid):
                    r1_type, r1_clip, r1_clip_type = CIGAR_field_classifier(r1_cigar_raw,experiment_type,contact[header_dict[r1_strand_header]])
                    r2_type, r2_clip, r2_clip_type = CIGAR_field_classifier(r2_cigar_raw,experiment_type,contact[header_dict[r2_strand_header]])
                    write_validation_reject(
                        CIGAR_filtered_out_contacts,
                        line,
                        mode,
                        [ r1_type, str(sum(r1_clip)), r1_clip_type, r2_type, str(sum(r2_clip)), r2_clip_type],
                        validation_reject_statistics,
                        "invalid_NM"
                    )
                    continue

                # Classify CIGARs and calculate final edit distances
                r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, \
                r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = \
                    CIGAR_field_classifier_parent(
                        pairtype,
                        edit_dist_type,
                        r1_cigar_raw,
                        contact[header_dict[r1_strand_header]],
                        r1_NM_raw,
                        r2_cigar_raw,
                        contact[header_dict[r2_strand_header]],
                        r2_NM_raw,
                        experiment_type
                    )
                #Apply experiment-specific splice-gap rules
                if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
                    r1_N_ok = validate_cigar_N(r1_cigar_raw,experiment_type,"rna")
                    r2_N_ok = validate_cigar_N(r2_cigar_raw,experiment_type,"dna")
                elif experiment_type == 'OTA_PE':
                    r1_N_ok = validate_cigar_N(r1_cigar_raw,experiment_type,"dna1")
                    r2_N_ok = validate_cigar_N(r2_cigar_raw,experiment_type,"dna2")
                else:
                    r1_N_ok = True
                    r2_N_ok = True

                # Validate MAPQ values
                r1_mapQ_raw = contact[header_dict[r1_mapQ_header]]
                r2_mapQ_raw = contact[header_dict[r2_mapQ_header]]
                r1_mapQ_valid = validate_MAPQ(r1_mapQ_raw)
                r2_mapQ_valid = validate_MAPQ(r2_mapQ_raw)
                if not (r1_mapQ_valid and r2_mapQ_valid):
                    write_validation_reject(
                        CIGAR_filtered_out_contacts,
                        line,
                        mode,
                        cigar_explorer_fields( r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type ),
                        validation_reject_statistics,
                        "invalid_MAPQ"
                    )
                    continue

                r1_mapQ, r2_mapQ = (int(r1_mapQ_raw), int(r2_mapQ_raw))

                # Validate genomic coordinates
                r1_start_raw = contact[header_dict[r1_start_header]]
                r1_end_raw = contact[header_dict[r1_end_header]]
                r2_start_raw = contact[header_dict[r2_start_header]]
                r2_end_raw = contact[header_dict[r2_end_header]]

                r1_coordinates_valid = validate_coordinates(r1_start_raw, r1_end_raw)
                r2_coordinates_valid = validate_coordinates(r2_start_raw, r2_end_raw)
                if not (r1_coordinates_valid and r2_coordinates_valid):
                    write_validation_reject(
                        CIGAR_filtered_out_contacts,
                        line,
                        mode,
                        cigar_explorer_fields( r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type ),
                        validation_reject_statistics,
                        "invalid_coordinates"
                    )
                    continue

                r1_start = int(r1_start_raw)
                r1_end = int(r1_end_raw)
                r2_start = int(r2_start_raw)
                r2_end = int(r2_end_raw)

                # Check paired-end genomic distance for OTA_PE
                distanceIsNormal = True if experiment_type in [
                    'ATA, not iMARGI',
                    'ATA, iMARGI'
                ] else check_distance_and_chromosomes(
                    contact[header_dict[r1_chr_header]],
                    r1_start,
                    r1_end,
                    contact[header_dict[r2_chr_header]],
                    r2_start,
                    r2_end,
                    OTA_PE_distanceThreshold
                )

                # Apply final filtering criteria
                if (
                    (r1_final_edit_dist <= r1_finalThreshold_edit_dist)
                    and (r2_final_edit_dist <= r2_finalThreshold_edit_dist)
                    and distanceIsNormal
                    and r1_N_ok
                    and r2_N_ok
                    and (r1_mapQ >= r1_mapQ_threshold)
                    and (r2_mapQ >= r2_mapQ_threshold)
                ):
                    output = output_string_modifier(line, mode, header_dict, edit_dist_type, experiment_type, N_softClipp_bp_or_NM_for_filtered_data, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type)
                    CIGAR_filtered_contacts.write(output + "\n")
                    if ucaRNA_enabled:
                        id_reads_for_ucaRNAs.write(contact[header_dict[SRR_ID_header]] + "\t" + pairtype + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)
                else:
                    output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, r2_cigar_type, str(sum(r2_N_softClipp_bp)), r2_softClipp_type])
                    CIGAR_filtered_out_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_out_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)

            else: # r1-only processing: ATA non-UU contacts, OTA_SE, RNAseq_SE, and RNAseq_PE
                # Validate r1 CIGAR before classification
                r1_cigar_raw = contact[header_dict[r1_cigar_header]]
                r1_cigar_valid = parse_cigar_tokens(r1_cigar_raw) is not None
                if not r1_cigar_valid:
                    write_validation_reject( CIGAR_filtered_out_contacts, line, mode, ["invalid","*","*","*","*","*"], validation_reject_statistics, "invalid_CIGAR_r1" )
                    continue

                # Validate NM before edit-distance calculation
                r1_NM_raw = contact[header_dict[r1_NM_header]]
                r1_NM_valid = validate_NM(r1_NM_raw)
                if not r1_NM_valid:
                    r1_type, r1_clip, r1_clip_type = CIGAR_field_classifier(r1_cigar_raw, experiment_type, contact[header_dict[r1_strand_header]])
                    write_validation_reject( CIGAR_filtered_out_contacts, line, mode, [r1_type,str(sum(r1_clip)),r1_clip_type,"*","*","*"], validation_reject_statistics, "invalid_NM" )
                    continue

                # Classify CIGAR and calculate final edit distance
                r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r1_final_edit_dist, \
                r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type, r2_final_edit_dist = \
                    CIGAR_field_classifier_parent(
                        pairtype,
                        edit_dist_type,
                        r1_cigar_raw,
                        contact[header_dict[r1_strand_header]],
                        r1_NM_raw,
                        contact[header_dict[r2_cigar_header]],
                        contact[header_dict[r2_strand_header]],
                        contact[header_dict[r2_NM_header]],
                        experiment_type
                    )

                # Apply experiment-specific splice-gap rules
                if experiment_type in ['ATA, not iMARGI', 'ATA, iMARGI']:
                    read_role = "rna"
                elif experiment_type == 'OTA_SE':
                    read_role = "dna"
                elif experiment_type == 'RNAseq_SE':
                    read_role = "rna"
                elif experiment_type == 'RNAseq_PE':
                    read_role = "rna1"
                else:
                    raise ValueError(f"Unsupported experiment type: {experiment_type}")

                r1_N_ok = validate_cigar_N(r1_cigar_raw,experiment_type,read_role)

                # Validate MAPQ
                r1_mapQ_raw = contact[header_dict[r1_mapQ_header]]
                r1_mapQ_valid = validate_MAPQ(r1_mapQ_raw)
                if not r1_mapQ_valid:
                    write_validation_reject(
                        CIGAR_filtered_out_contacts,
                        line,
                        mode,
                        cigar_explorer_fields(r1_cigar_type,r1_N_softClipp_bp,r1_softClipp_type),
                        validation_reject_statistics,
                        "invalid_MAPQ"
                    )
                    continue

                r1_mapQ = int(r1_mapQ_raw)

                # Validate genomic coordinates
                r1_coordinates_valid = validate_coordinates(contact[header_dict[r1_start_header]], contact[header_dict[r1_end_header]])
                if not r1_coordinates_valid:
                    write_validation_reject(
                        CIGAR_filtered_out_contacts,
                        line,
                        mode,
                        cigar_explorer_fields(r1_cigar_type,r1_N_softClipp_bp,r1_softClipp_type),
                        validation_reject_statistics,
                        "invalid_coordinates"
                    )
                    continue

                # Apply final filtering criteria
                if (
                    (r1_final_edit_dist <= r1_finalThreshold_edit_dist)
                    and r1_N_ok
                    and (r1_mapQ >= r1_mapQ_threshold)
                ):
                    output = output_string_modifier(line, mode, header_dict, edit_dist_type, experiment_type, N_softClipp_bp_or_NM_for_filtered_data, SRR_ID_header, pairtype_header, r1_chr_header, r1_start_header, r1_end_header, r1_strand_header, r1_cigar_header, r1_secondary_alignments_header, r1_other_tags_header, r1_NM_header, r2_chr_header, r2_start_header, r2_end_header, r2_strand_header, r2_cigar_header, r2_secondary_alignments_header, r2_other_tags_header, r2_NM_header, r1_mapQ_header, r2_mapQ_header, r1_cigar_type, r1_N_softClipp_bp, r1_softClipp_type, r2_cigar_type, r2_N_softClipp_bp, r2_softClipp_type)
                    CIGAR_filtered_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)
                    if ucaRNA_enabled:
                        id_reads_for_ucaRNAs.write(contact[header_dict[SRR_ID_header]] + "\t" + pairtype + "\n")
                else:
                    output = line if mode != "explorer" else line + "\t" + "\t".join([r1_cigar_type, str(sum(r1_N_softClipp_bp)), r1_softClipp_type, "*", "*", "*"])
                    CIGAR_filtered_out_contacts.write(output + "\n")
                    CIGAR_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_cigar_type, r2_cigar_type)
                    softClipp_NM_statistics(N_softClipp_bp_or_NM_for_filtered_out_data, r1_N_softClipp_bp, r2_N_softClipp_bp, contact[header_dict[r1_NM_header]], contact[header_dict[r2_NM_header]], r1_chr_header, r2_chr_header)


    # Closing files
    raw_contacts.close()
    CIGAR_filtered_contacts.close()
    CIGAR_filtered_out_contacts.close()
    if ucaRNA_enabled:
        id_reads_for_ucaRNAs.close()

    save_cigar_statistics(CIGAR_statistics_for_filtered_data, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, "filtered")
    save_cigar_statistics(CIGAR_statistics_for_filtered_out_data, experiment_type, r1_chr_header, r2_chr_header, output_file_path, input_file_name, "out")
    save_validation_reject_statistics(validation_reject_statistics, output_file_path, input_file_name)

    plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_for_filtered_data, experiment_type, r1, r2, 'filtered_data', output_file_path)
    plot_N_softClipp_or_NM(N_softClipp_bp_or_NM_for_filtered_out_data, experiment_type, r1, r2, 'out_data', output_file_path)



# edit_dist_type              = 'NM + N_softClipp_bp' or 'NM'
# r1_finalThreshold_edit_dist = integer >= 0
# r2_finalThreshold_edit_dist = integer >= 0
# r1_mapQ_threshold           = integer >= 0
# r2_mapQ_threshold           = integer >= 0
# OTA_PE_distanceThreshold    = integer >= 0
# Assembly_of_ucaRNAs         = 'yes' or 'no'
# mode                        = 'not explorer' or 'explorer'
# experiment_type             = 'ATA, not iMARGI', 'ATA, iMARGI', 'OTA_SE', 'OTA_PE', 'RNAseq_SE', 'RNAseq_PE'
# input_file_name             = input contacts-table filename
# input_file_path             = directory containing the input file
# output_file_path            = output directory

#/usr/bin/time -v python EditDistance_CIGAR_filter.py "NM + N_softClipp_bp" 1 1 0 0 200 "yes" "not explorer" "ATA, not iMARGI" "SRR17331267_K562_CTCF_RedChIP_rep1_Unique_RNA.tab.rc" "/path/to/input/file/" "/path/to/output/dir/"
# editDistance_and_CIGAR_filter(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), sys.argv[7], sys.argv[8], sys.argv[9], sys.argv[10], sys.argv[11], sys.argv[12])
if __name__ == "__main__":

    if len(sys.argv) != 13:
        raise SystemExit(
            "Usage:\n"
            "python EditDistance_CIGAR_filter.py "
            "<edit_dist_type> "
            "<r1_edit_distance_threshold> "
            "<r2_edit_distance_threshold> "
            "<r1_MAPQ_threshold> "
            "<r2_MAPQ_threshold> "
            "<OTA_PE_distance_threshold> "
            "<Assembly_of_ucaRNAs> "
            "<mode> "
            "<experiment_type> "
            "<input_file_name> "
            "<input_file_path> "
            "<output_file_path>"
        )

    edit_dist_type = sys.argv[1]
    r1_finalThreshold_edit_dist = int(sys.argv[2])
    r2_finalThreshold_edit_dist = int(sys.argv[3])
    r1_mapQ_threshold = int(sys.argv[4])
    r2_mapQ_threshold = int(sys.argv[5])
    OTA_PE_distanceThreshold = int(sys.argv[6])
    Assembly_of_ucaRNAs = sys.argv[7]
    mode = sys.argv[8]
    experiment_type = sys.argv[9]
    input_file_name = sys.argv[10]
    input_file_path = sys.argv[11]
    output_file_path = sys.argv[12]

    if edit_dist_type not in ["NM", "NM + N_softClipp_bp"]:
        raise ValueError(
            f"Unsupported edit distance type: {edit_dist_type}"
        )

    if Assembly_of_ucaRNAs not in ["yes", "no"]:
        raise ValueError(
            "Assembly_of_ucaRNAs must be 'yes' or 'no'"
        )

    if mode not in ["not explorer", "explorer"]:
        raise ValueError(
            "mode must be 'not explorer' or 'explorer'"
        )

    valid_experiment_types = [
        "ATA, not iMARGI",
        "ATA, iMARGI",
        "OTA_SE",
        "OTA_PE",
        "RNAseq_SE",
        "RNAseq_PE"
    ]

    if experiment_type not in valid_experiment_types:
        raise ValueError(
            f"Unsupported experiment type: {experiment_type}"
        )

    if any(
        value < 0
        for value in [
            r1_finalThreshold_edit_dist,
            r2_finalThreshold_edit_dist,
            r1_mapQ_threshold,
            r2_mapQ_threshold,
            OTA_PE_distanceThreshold
        ]
    ):
        raise ValueError(
            "Threshold values must be non-negative integers"
        )

    editDistance_and_CIGAR_filter(
        edit_dist_type,
        r1_finalThreshold_edit_dist,
        r2_finalThreshold_edit_dist,
        r1_mapQ_threshold,
        r2_mapQ_threshold,
        OTA_PE_distanceThreshold,
        Assembly_of_ucaRNAs,
        mode,
        experiment_type,
        input_file_name,
        input_file_path,
        output_file_path
    )