#!/bin/bash

set -ex

if [ ! -f 'bmark.ok' ]; then

    $(git rev-parse --show-toplevel)/benchmarking/bmark_nb_runner.py \
        --analysisType ref_guided --truth_gtf reference_data/minigenome.UTRs_trimmed_1isoformref_expressed.gtf --truth_quant reference_data/morf2_ont_merged_annot_compat_sorted_tn_counts.tsv --truth_reduced_gtf reference_data/minigenome.UTRs_trimmed_1isoformref_reduced.gtf

    touch bmark.ok

fi
