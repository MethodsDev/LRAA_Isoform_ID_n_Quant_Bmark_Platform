#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType quant_only --truth_gtf ../../reference_data/SIRV_isoforms_multi-fasta-annotation_C_170612a.gtf --truth_quant ../../reference_data/BT474_E0_merged_sirv_sorted_groundtruth_E0.tsv

    touch bmark.ok
    
fi
