#!/bin/bash

set -e

tests=(
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
)

for test in "${tests[@]}"
do
    echo
    echo "============================================================"
    echo "RUNNING: $test"
    echo "============================================================"

    python "$test"
done

echo
echo "============================================================"
echo "ALL TEST SCRIPTS COMPLETED"
echo "============================================================"