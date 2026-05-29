#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType ref_guided --truth_gtf ../../reference_data/SIRV_isoforms_multi-fasta-annotation.expressed.gtf --truth_reduced_gtf ../../reference_data/SIRV_isoforms_multi-fasta-annotation.reduced_1isoform_pigeon.gtf --truth_quant ../../reference_data/UHRR_E0_merged_sirv_sorted_groundtruth_E0.tsv

    touch bmark.ok

fi

